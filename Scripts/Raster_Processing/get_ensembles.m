% Function to Indentify Ensembles in Neural Activity
% by Clustering (hierarchichally) Neural *Coactivity* 
% Applies to the Raster after Selection 
% ***   OUTPUT has to be named with the sufix: _Analysis *****
% Neither need of tangled abstract stuff nor Monte Carlos Simulations
% Input
%   R:          Raster of Neural Activity (Matrix Cells x Frames)
% Output
%   R_Analysis:  Strcuture of Neural Ensemble Analysis
function R_Analysis = get_ensembles(R)
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
    R_Analysis.Data.Data=R'; % Frames x Cells (compatibility wit JP's NN)
    fprintf('>> Finishing script...');
    Analyze=false;
else
    % Just Initialize:
    R_Analysis.Data.Data=R';  % Frames x Cells (compatibility wit JP's NN)
    R_Analysis.Peaks.Threshold=0;               %CAG_TH
    R_Analysis.Clustering.TotalStates=0;        % Nensemebles
    R_Analysis.Peaks.Index=zeros(1,Frames);     % Set 1 to Frames(CAG=>Threshold)
    R_Analysis.Clustering.VectorStateIndex=[];  % Sequence of Ensembles: frame_ensembles
    Analyze=true;
end
if Analyze
    fprintf('> Maximum Coactive Neurons   %i  \n',MaxCAG);
    % CAGThresholdValues=1:MaxCAG;
    % About the CLustering Method
    SimMethod='hamming'; % From Binary Codes used in Communications
    % This will increase as ensembles ain't a complete subset from another:
    disp('> Analysis Started...')
    %% MAIN LOOPs
    % EnsembleOK='true';
    % Nensembles=2;
    NensemblesTotal=10;
    Dunn=zeros(MaxCAG,NensemblesTotal);
    PercAN=zeros(1,MaxCAG);
    for Nensembles=2:NensemblesTotal
        for CAGindex=1:MaxCAG
            % Clustering magic
            NeuroVectors=zeros(numel(ActiveCells),Nensembles);
            EnsembledNeurons={};
            fprintf('>>> Clustering for %i Ensembles ',Nensembles);
            fprintf('& %i Coactive Neurons\n',CAGindex);
            Rclust=Ractive(:,CAG>=CAGindex);
            if Nensembles<3
                % Compute just ONCE
                ActiveCellsClust=find(sum(Rclust,2));
                PercAN(CAGindex)=numel(ActiveCellsClust)/numel(ActiveCells);
                fprintf('Active Neurons: %3.3f \n',PercAN(CAGindex));
            end 
            [~,ActiveFrames]=size(Rclust);
            if ActiveFrames>Nensembles
                Distance=squareform(pdist(Rclust',SimMethod));
                Sim=1-Distance;
                LinkageMethod=HBestTree_JPplus(Sim);    % Output
                Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
                frame_ensembles=cluster(Tree,'maxclust',Nensembles); % Output
                % Get Neurons @ Ensembles:
                for n=1:Nensembles
                    EnsembledNeurons{n}=find(sum(Rclust(:,frame_ensembles==n),2));
                    NeuroVectors(EnsembledNeurons{n},n)=1;
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
    disp('>> Scanning Clusters: Completed.')
    %% CAG Threshold *******************************
    ActiveNeuronsRation=1;
    CAGwithAN=find(PercAN>=ActiveNeuronsRation);
    CAGthreshold=CAGwithAN(end);
    PercAN(CAGthreshold:MaxCAG)
    %% PLotting Result
    disp('>> Plotting Results:')
    ScanClusterFig=figure('name',['Scanning Ensembles for ',inputname(1)],'numbertitle','off');
    ScanClusterFig.Position=[778 264 576 399];
    H1ax=subplot(3,1,[1,2]);
    imagesc(Dunn); colormap(gca,jet);
    ScanClusterFig.Children.YDir='normal';
    ScanClusterFig.Children.YTick=1:MaxCAG;
    colorbar('peer',gca);
    H2ax=subplot(3,1,3);
    H2ax.XLim=[0.5,NensemblesTotal+0.5];
    H2ax.Position=[0.13,0.14,0.66,0.19];
    bar(H2ax,1:NensemblesTotal,max(Dunn),'k')
    xlabel(H2ax,'Number of Ensembles')
    ylabel(H1ax,'CAG Threshold')
    ylabel(H2ax,'maxDunn')
    title(H1ax,'Index Dunn')
    %% Find the BEST Clusterization
    MaxDunn=max(max(Dunn));
    [rowM,colM]=find(Dunn==MaxDunn);
    if numel(rowM)==1
        CAG_TH=rowM;
        Nensembles=colM;
        fprintf('>>Global Max Dunn Index Found\n>>With %i Ensembles & for %i Coactive Neurons\n',Nensembles,CAG_TH)
    else
        disp('>>Local Maximum Dunn Index Found');
        disp('>>Searching minimum Number of Ensembles and Maximum CAG Threshold...')
        Nensembles=min(colM);
        CAG_TH=max(rowM(colM==Nensembles));
    end
    fprintf('>> Clustering with %i Ensembles & for %i Coactive Neurons\n',Nensembles,CAG_TH)
    Rclust=Ractive(:,CAG>=CAG_TH);
    Distance=squareform(pdist(Rclust',SimMethod));
    Sim=1-Distance;
    LinkageMethod=HBestTree_JPplus(Sim);    % Output
    Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
    frame_ensembles=cluster(Tree,'maxclust',Nensembles); % Output
    % OUTPUT
    R_Analysis.Data.Data=R';  % Frames x Cells (compatibility wit JP's NN)
    R_Analysis.Peaks.Threshold=CAG_TH;
    R_Analysis.Clustering.TotalStates=Nensembles;
    R_Analysis.Clustering.Tree=Tree;
    R_Analysis.Clustering.Linkage=LinkageMethod;
    R_Analysis.Peaks.Index(CAG>=CAG_TH)=1;
    R_Analysis.Clustering.VectorStateIndex=frame_ensembles;  % Sequence of Ensembles: frame_ensembles
    % PLOt cluster
    hold (H1ax,'on'); hold (H2ax,'on');
    plot(H1ax,Nensembles,CAG_TH,'Marker','pentagram','Color','g','MarkerFaceColor','g','MarkerSize',13);
    plot(H2ax,Nensembles,Dunn(Nensembles,CAG_TH),'Marker','pentagram','Color','g','MarkerFaceColor','g','MarkerSize',10);
    H2ax.XLim=[0.5,NensemblesTotal+0.5];
    hold (H1ax,'off'); hold (H2ax,'off');
end
fprintf('\n>>Script to search Neural Ensembles has ended.\n')