% Function to Indentify Ensembles in Neural Activity
% by Clustering (hierarchichally) Neural *Coactivity* 
% Applies to the Raster after Selection 
% ***   OUTPUT has to be named with the sufix: _Analysis *****
% Neither need of tangled abstract stuff nor Monte Carlos Simulations
% Input
%   R:          Raster of Neural Activity (Matrix Cells x Frames)
% Output
%   R_Analysis: Raster Structur with Fields:
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
function R_Analysis = get_bayes_ensembles(R)
%% Setup
% About the Raster
disp('> Getting data...')

[Cells,Frames]=size(R); % Always More Frames than Cells
if Cells>Frames
    R=R';
    [Cells,Frames]=size(R);
    disp('>>Raster Transposed.')
end
ActiveCells=find(sum(R,2));
% GET ONLY ACTIVE CELLS!
Ractive=R(ActiveCells,:);
CAG=sum(Ractive);
MaxCAG=max(CAG);
%% AS IF ITS EMPTY *****??????
% Initialize Output
if MaxCAG==0
    fprintf('\n>> Raster with Zero Activity\n\n');
    R_Analysis.Data.Data=R; % Frames x Cells (compatibility wit JP's NN)
    fprintf('>> Finishing script...');
    Analyze=false;
else
    % Just Initialize:
    R_Analysis.Data.Data=R';  % Frames x Cells (compatibility wit JP's NN)
    R_Analysis.Peaks.Threshold=0;               %CAG_TH
    R_Analysis.Clustering.TotalStates=0;        % Nensemebles
    R_Analysis.Peaks.Index=zeros(1,Frames);     % Set 1 to Frames(CAG=>Threshold)
    R_Analysis.Clustering.VectorStateIndex=[];  % Sequence of Ensembles: frame_ensembles
    R_Analysis.Classifier.Model=[];             % Naive Bayes Classifier
    R_Analysis.Classifier.ValidationError=[];   % Validation Error
    Analyze=true;
end
%% START ANALYSIS
if Analyze
    fprintf('> Maximum Coactive Neurons   %i  \n',MaxCAG);
    % CAGThresholdValues=1:MaxCAG;
    % About the CLustering Method
    SimMethod='hamming'; % From Binary Codes used in Communications
    % This will increase as ensembles ain't a complete subset from another:
    disp('> Analysis Initialized ...')
    %% GLOBAL ANALYSIS SETUP
    tic;
    % GET Possible CAG Thresholds:
    ActiveNeuronsRatio=0.75;    % INPUT
    ActiveTime=0.5;             % INPUT
    CAGwithAN=[];
    for CAGindex=1:MaxCAG
        Rclust=Ractive(:,CAG>=CAGindex);
        ActiveCellsClust=find(sum(Rclust,2));
        PercAN(CAGindex)=numel(ActiveCellsClust)/numel(ActiveCells);
        ActTime(CAGindex)=sum((sum(Rclust)>0))/size(Ractive,2);
        fprintf('For %i CA Neurons-> %3.1f%% Active Neurons %3.1f%% of the Time\n',CAGindex,100*PercAN(CAGindex),100*ActTime(CAGindex));
        if PercAN(CAGindex)>=ActiveNeuronsRatio && ActTime(CAGindex)>=ActiveTime
            CAGwithAN=[CAGwithAN,CAGindex];
        end
    end
    % Corrections --------------------
    if numel(CAGwithAN)==1
        if CAGwithAN==1 && CAGwithAN<MaxCAG
            CAGwithAN=[CAGwithAN,2];
        elseif CAGwithAN<MaxCAG
            CAGwithAN=[CAGwithAN,CAGwithAN+1];
        end
    end
    % If there is not too much active time, use
    % Only Ratio of Active Neurons
    if isempty(CAGwithAN)    
        CAGwithAN=find(PercAN>ActiveNeuronsRatio);
    end
    % If still EMPTY->>Accept ALL
    if isempty(CAGwithAN)
        CAGwithAN=1:MaxCAG;
    end
    
    % Ensembles Setup ************
    NensemblesTotal=5;
    Dunn=zeros(numel(CAGwithAN),NensemblesTotal);
    ClusterOK=zeros(numel(CAGwithAN),NensemblesTotal);
    MinIntraClusterSim=zeros(numel(CAGwithAN),NensemblesTotal);
    % ClassRatio=zeros(numel(CAGwithAN),NensemblesTotal);
    ErrorClass=ones(numel(CAGwithAN),NensemblesTotal);
    %% MAIN LOOPS
    for Nensembles=2:NensemblesTotal
        for CAGindex=1:numel(CAGwithAN)
            % Clustering magic
            NeuroVectors=zeros(numel(ActiveCells),Nensembles);
            EnsembledNeurons={};
            fprintf('>>> Clustering for %i Ensembles ',Nensembles);
            fprintf('& %i Coactive Neurons\n',CAGwithAN(CAGindex));
            Rclust=Ractive(:,CAG>=CAGwithAN(CAGindex));
            [~,ActiveFrames]=size(Rclust);
            if ActiveFrames>Nensembles
                Distance=squareform(pdist(Rclust',SimMethod));
                Sim=1-Distance;
                LinkageMethod=HBestTree_JPplus(Sim);    % Output
                Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
                frame_ensembles=cluster(Tree,'maxclust',Nensembles); % Output
                % Evaluate w/Naive Bayes Classiffier:
                [~,ECV]=Nbayes_Ensembles(Rclust,frame_ensembles);
                ErrorClass(CAGindex,Nensembles)=ECV;
                % Get Neurons @ Ensembles:
                for n=1:Nensembles
                    EnsembledNeurons{n}=find(sum(Rclust(:,frame_ensembles==n),2));
                    NeuroVectors(EnsembledNeurons{n},n)=1;
                    % Intra Cluster Distances
                    Rclustens=Rclust(EnsembledNeurons{n},frame_ensembles==n);
                    ClustDist=1-pdist(Rclustens',SimMethod);
                    if isempty(ClustDist(ClustDist>0))
                        MinCLustSimm(n)=0;
                    else
                        MinCLustSimm(n)=min(ClustDist(ClustDist>0));
                    end
                end
                MinIntraClusterSim(CAGindex,Nensembles)=min(MinCLustSimm);
                % Criteria to exclude Clusterings
                % A ensembles belongs to other
                EnsemblesIndx=isinensemble(NeuroVectors);
                if isempty(EnsemblesIndx)
                    % ACCEPTED CLUSTERING
                    ClusterOK(CAGindex,Nensembles)=1;
                else
                    % REJECTED CLUSTERING
                    disp('>>Redundant Ensemble(s) Found!!!')
                    [NRepEns,~]=size(EnsemblesIndx);
                    for nE=1:NRepEns
                        fprintf('\nEnsemble: %i is in Ensemble %i\n',EnsemblesIndx(nE,1),EnsemblesIndx(nE,2));
                    end
                end
                % 1-neuron Ensembles
                if sum(sum(NeuroVectors)==1)>0
                    disp('1-neuron Ensembles');
                    ClusterOK(CAGindex,Nensembles)=0;
                end
                % DunnIndex:
                Dhamm=pdist(NeuroVectors',SimMethod);
                if isempty(Dhamm); Dhamm=0; end;
                Lens=sum(NeuroVectors)/numel(ActiveCells);
                if isempty(Lens)
                    Dunn(CAGindex,Nensembles)=0;
                else
                    Dunn(CAGindex,Nensembles)=min(Dhamm)/max(Lens); % min distance among ensembles divided by maximum length of ensembles
                end
            else
                disp('>>>> Not enough frames to Analyze.')
                Dunn(CAGindex,Nensembles)=0;
            end
        end
        fprintf('\n\n');
    end
    toc;
    %% GET BEST CLUSTERING:
    % Minimum DIstance INTRA cluster
    % Best Classification
    % & Clusterin Index OK->true
%     %% OPTION A
%     IntraClusterDistances=unique(sort(MinIntraClusterSim(:),'descend'),'stable');
%     CheckCluster=true; n=1;
%     while CheckCluster
%         [CAGs,Enss]=find(MinIntraClusterSim==IntraClusterDistances(n));
%         % Get Classification Errors & ClusterOK indicator
%         ClassError=[];
%         okclust=[];
%         for c=1:numel(CAGs)
%             ClassError(c)=ErrorClass(CAGs(c),Enss(c));
%             okclust(c)=ClusterOK(CAGs(c),Enss(c));
%         end
%         [~,Okindx]=min(ClassError(okclust>0));
%         if isempty(Okindx)
%             disp('>>Keep Searching...>>')
%             n=n+1;
%             % Getting Next "Best Clusters"
%         else
%             disp('>>Ready')
%             CAG_THindx=CAGs(Okindx);
%             Nensembles=Enss(Okindx);
%             CheckCluster=false;
%         end
%     end
%     CAG_TH=CAGwithAN(CAG_THindx);
    %% OPTION B
%     ErroClassAll=sort(ErrorClass(1,:)); % Taking into account all the Raster
%     CheckCluster=true; n=1;
%     while CheckCluster
%         % Ensmables that CLassify best ALL the activity
%         Enss=find(ErrorClass(1,:)==ErroClassAll(n),1);
%         % Get best Clusterin CAG threshold
%         [~,CAGs]=max(MinIntraClusterSim(:,Enss));
%         
%         if ClusterOK(CAGs,Enss)==0
%             disp('>>Keep Searching...>>')
%             n=n+1;
%             % Getting Next "Best Clusters"
%         else
%             disp('>>Ready')
%             CAG_THindx=CAGs;
%             Nensembles=Enss;
%             CheckCluster=false;
%         end
%     end
%     CAG_TH=CAGwithAN(CAG_THindx);
    
    %% OPTION C
    ErroClassAll=sort(ErrorClass(:)); % Taking into account all the Raster
    CheckCluster=true; n=1;
    while CheckCluster
        % Ensmables that CLassify best ALL the activity
        [CAGs,Enss]=find(ErrorClass==ErroClassAll(n),1);
        % Get best Clusterin CAG threshold
        % [~,CAGs]=max(MinIntraClusterSim(:,Enss));
        
        if ClusterOK(CAGs,Enss)==0
            disp('>>Keep Searching...>>')
            n=n+1;
            % Getting Next "Best Clusters"
        else
            disp('>>Ready')
            CAG_THindx=CAGs;
            Nensembles=Enss;
            CheckCluster=false;
        end
        if n>numel(ErroClassAll)
            CheckCluster=false;
        end
    end
    CAG_TH=CAGwithAN(CAG_THindx);
    
    % Get CAG Threshold: The one that uses is not affected by
    % number of ensembles
    % disp('Getting best Classification:')
    %[CAG_THindx,Nensembles]=find(ErrorClass==min(ErrorClass(ClusterOK>0)));
    
    %MeanErrorNensembles=mean(ErrorClass,2);
    %[~,CAG_THindx]=min(MeanErrorNensembles(MeanErrorNensembles>0));
    % Get Number of Ensembles: 
    % [~,Nensembles]=findpeaks(-ErrorClass(CAG_TH,:),'Npeaks',1);
    % if isempty(Nensembles)
    %    disp('>>>Atypical Clustering')
    %    [~,Nensembles]=findpeaks(-abs(ErrorClass(CAG_TH,:)-mean(ErrorClass(CAG_TH,2:end))),'Npeaks',1);
    % end
    
    %% Plotting Result
    disp('>> Plotting Results:')
    ScanClusterFig=figure('name',['Clustering Analysis for ',inputname(1)],'numbertitle','off');
    ScanClusterFig.Position=[778 264 576 399];
    H1ax=subplot(3,1,[1,2,3]);
    imagesc(ErrorClass); colormap(gca,jet);
    ScanClusterFig.Children.YDir='normal';
    ScanClusterFig.Children.YTick=1:MaxCAG;
    colorbar('peer',gca);
    % H2ax=subplot(3,1,3);
    % H2ax.XLim=[0.5,NensemblesTotal+0.5];
    % H2ax.Position=[0.13,0.14,0.66,0.19];
    % bar(H2ax,1:NensemblesTotal,1,'k')
    xlabel(H1ax,'Number of Ensembles')
    ylabel(H1ax,'CAG Threshold')
    % ylabel(H2ax,'Classification Error')
    title(H1ax,'Validation Error')
    hold (H1ax,'on'); 
    plot(H1ax,Nensembles,CAG_TH,'Marker','pentagram','Color','g','MarkerFaceColor','g','MarkerSize',13);
    hold (H1ax,'off');
    %% Find the BEST Clusterization
%     MaxDunn=max(max(Dunn));
%     [rowM,colM]=find(Dunn==MaxDunn);
%     if numel(rowM)==1
%         CAG_TH=rowM;
%         Nensembles=colM;
%         fprintf('>>Global Max Dunn Index Found\n>>With %i Ensembles & for %i Coactive Neurons\n',Nensembles,CAG_TH)
%     else
%         disp('>>Local Maximum Dunn Index Found');
%         disp('>>Searching minimum Number of Ensembles and Maximum CAG Threshold...')
%         Nensembles=min(colM);
%         CAG_TH=max(rowM(colM==Nensembles));
%     end
    %% Final Clustering
    fprintf('>> Clustering with %i Ensembles & for %i Coactive Neurons\n',Nensembles,CAG_TH)
    % CLuster and plot SimmMat&dendro:
    R_Analysis=raster_cluster(R,CAG_TH,Nensembles,SimMethod,1);
%     Rclust=Ractive(:,CAG>=CAG_TH);
%     Distance=squareform(pdist(Rclust',SimMethod));
%     Sim=1-Distance;
%     LinkageMethod=HBestTree_JPplus(Sim);    % Output
%     Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
%     frame_ensembles=cluster(Tree,'maxclust',Nensembles); % Output
%     % RE-LABEL Ensembles
%     AppearSequence=unique(frame_ensembles,'stable');
%     relabel_frame_ensembles=zeros(size(frame_ensembles));
%     for n=1:Nensembles
%         relabel_frame_ensembles(frame_ensembles==AppearSequence(n))=n;
%     end
%     % Evaluate w/Naive Bayes Classiffier:
%     [Model,ECV]=Nbayes_Ensembles(Rclust,relabel_frame_ensembles);
%     % OUTPUT
%     R_Analysis.Data.Data=R';  % Frames x Cells (compatibility wit JP's NN)
%     R_Analysis.Peaks.Threshold=CAG_TH;
%     R_Analysis.Clustering.TotalStates=Nensembles;
%     R_Analysis.Clustering.Tree=Tree;
%     R_Analysis.Clustering.Linkage=LinkageMethod;
%     R_Analysis.Peaks.Index(CAG>=CAG_TH)=1;
%     R_Analysis.Clustering.VectorStateIndex=relabel_frame_ensembles;  % Sequence of Ensembles: frame_ensembles
%     R_Analysis.Classifier.Model=Model;             % Naive Bayes Classifier
%     R_Analysis.Classifier.ValidationError=ECV;      % Validation Error
end
fprintf('>> Clustering with %i Ensembles \n& for %i Coactive Neurons\n',Nensembles,CAG_TH)
fprintf('\n>>Script to search Neural Ensembles has ended.\n')