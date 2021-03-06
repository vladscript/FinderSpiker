% Hierachichal Clustering of Raster Coactivity
% .
% Input
%   R:          Raster of Neural Activity(Frames x Cells Raster)
%   CAG_TH:     Coactivity Threshold
%   Nensembles: Number of Ensembles to Cluster
%   SimMethod:  Method of Simmilarity
% Ouput
%   Raster Structur with Fields:
%       Data.Data=R;  % Frames x Cells (compatibility wit JP's NN)
%       Peaks.Threshold=CAG_TH;
%       Clustering.TotalStates=Nensembles;
%       Clustering.Tree=Tree;
%       Clustering.Linkage=LinkageMethod;
%       Peaks.Index=zeros(1,numel(CAG));
%       Peaks.Index(CAG>=CAG_TH)=1;
%       Clustering.VectorStateIndex=frame_ensembles;  % Sequence of Ensembles: frame_ensembles
%       Classifier.Model=Model;             % Naive Bayes Classifier
%       Classifier.ValidationError=ECV;      % Validation Error

function R_Analysis = raster_cluster(R,CAG_TH,Nensembles,SimMethod,varargin)
if isempty(varargin)
    plotthings=false;
else
    plotthings=true;
end
[Frames,Cells]=size(R); % Always More Frames than Cells
if Cells>Frames
    R=R';
    % [Cells,Frames]=size(R);
    disp('>>Raster Transposed.')
end
ActiveCells=find(sum(R));
% GET ONLY ACTIVE CELLS!
Ractive=R(:,ActiveCells);
CAG=sum(Ractive,2);

MaxCAG=max(CAG);

if MaxCAG<CAG_TH
    fprintf('\n>> Raster with LOW Co-Activity \n\n');
    R_Analysis.Data.Data=R; % Frames x Cells (compatibility wit JP's NN)
    fprintf('>> Finishing script...');
    Analyze=false;
else
    Analyze=true;
end
if Analyze
    Rclust=Ractive(CAG>=CAG_TH,:);
    Distance=squareform(pdist(Rclust,SimMethod));
    Sim=1-Distance;
    LinkageMethod=HBestTree_JPplus(Sim);    % Output
    Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
    frame_ensembles=cluster(Tree,'maxclust',Nensembles); % Output
    % RE-LABEL Ensembles
    AppearSequence=unique(frame_ensembles,'stable');
    relabel_frame_ensembles=zeros(size(frame_ensembles));
    for n=1:Nensembles
        relabel_frame_ensembles(frame_ensembles==AppearSequence(n))=n;
    end
    % Plot Simmilarity Matrix & Dendrogram
    if plotthings
        figure;
        % [~,sortedframes]=sort(relabel_frame_ensembles);
        % Distance=squareform(pdist(Rclust(sortedframes,:),SimMethod));
        % Sim=1-Distance;
        Hmat=subplot(1,2,1);
        imagesc(Sim);
        Hden=subplot(1,2,2);
        dendrogram(Tree,0,'Orientation','right');
        Hmat.XTickLabel=[];
        % Hmat.YTickLabel=[];
        Hden.YTickLabel=[];
        Hden.XLabel.String=[LinkageMethod,' linkage'];
    end
    % Evaluate w/Naive Bayes Classiffier:
    [Model,ECV]=Nbayes_Ensembles(Rclust',relabel_frame_ensembles);
    % OUTPUT
    R_Analysis.Data.Data=R;  % Frames x Cells (compatibility wit JP's NN)
    R_Analysis.Peaks.Threshold=CAG_TH;
    R_Analysis.Clustering.TotalStates=Nensembles;
    R_Analysis.Clustering.Tree=Tree;
    R_Analysis.Clustering.Linkage=LinkageMethod;
    R_Analysis.Peaks.Index(CAG>=CAG_TH)=1;
    R_Analysis.Clustering.VectorStateIndex=relabel_frame_ensembles;  % Sequence of Ensembles: frame_ensembles
    R_Analysis.Classifier.Model=Model;             % Naive Bayes Classifier
    R_Analysis.Classifier.ValidationError=ECV;   
end
end
