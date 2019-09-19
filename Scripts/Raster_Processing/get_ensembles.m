% Function to Indentify Ensembles in Neural Activity
% by Clustering (hierarchichally) Neural *Coactivity* 
% 2 or more cells in experiment
% Ignores cells that fires less than 0.25% of the ensemble time
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
%       Clustering.Tree=Tree;   Hierarchichal Clustering
%       Clustering.Linkage=LinkageMethod;   Linkage
%       Peaks.Index=zeros(1,numel(CAG));
%       Peaks.Index(CAG>=CAG_TH)=1;
%       Clustering.VectorStateIndex=frame_ensembles;  % Sequence of Ensembles: frame_ensembles
%       Classifier.Model=Model;             % Naive Bayes Classifier
%       Classifier.ValidationError=ECV;      % Validation Error
function R_Analysis = get_ensembles(R,ThCAG,Nensambles)
%% Setup 
Load_Default_Clustering;

% About the Raster
disp('> Getting data ...')
[~,Frames]=size(R); % Always More Frames than Cells
ActiveCells=find(sum(R,2));
Ractive=R(ActiveCells,:);
CAG=sum(Ractive);
MaxCAG=max(CAG);
%% ASK IF ITS EMPTY *****??????
% Initialize Output
if MaxCAG==0 || MaxCAG<ThCAG
    fprintf('\n>> Raster with Zero or Low Activity\n\n');
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
    fprintf('> Maximum Coactive Neurons   %i  \n',MaxCAG,Nensambles);
    % CAGThresholdValues=1:MaxCAG;
    % This will increase as ensembles ain't a complete subset from another:
    disp('> Analysis Initialized ...')
    %% GLOBAL ANALYSIS SETUP
    tic;
    %     
    CAGwithAN=[];
    for CAGindex=2:MaxCAG
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
    % ClassRatio=zeros(numel(CAGwithAN),NensemblesTotal);
%     CAGwithAN=CAGwithAN(CAGwithAN>=ThCAG);
    ErrorClass=ones(numel(CAGwithAN),1);
    NensemblesOK=ones(numel(CAGwithAN),1);
    % CAGwithAN=CAGwithAN(CAGwithAN>=ThCAG);
    CAGwithAN=ThCAG;
    %% MAIN LOOPS CAG THReSHOLDS **************************************
    % ModelS=cell(numel(CAGwithAN),1);
    for CAGindex=1:numel(CAGwithAN)
        fprintf('>>> Clustering for  %i Coactive Neurons\n',CAGwithAN(CAGindex));
        OrgigiFrames{CAGindex}=find(CAG>=CAGwithAN(CAGindex));
        Rclust=Ractive(:,OrgigiFrames{CAGindex});
        [~,ActiveFrames]=size(Rclust);
        if ActiveFrames>1
            
            Distance=squareform(pdist(Rclust',SimMethod));
            Sim=1-Distance;
            LinkageMethod=HBestTree_JPplus(Sim);    % Output
            Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
            frame_ensembles=cluster(Tree,'maxclust',Nensambles); % START
            
            % [~,MaxDistanceIV]=nonsimilarframes(Rclust,SimMethod,0.5);
            % Check Ignored:
            signif_frames=1:size(Rclust(:,frame_ensembles>0),2);
            for t=1:3
                [~,ECV(t)]=Nbayes_Ensembles(Rclust(:,signif_frames),...
                    frame_ensembles(frame_ensembles>0));
            end
            % MaxDIV(CAGindex)=MaxDistanceIV;
            NensemblesOK(CAGindex)=numel(unique(frame_ensembles(frame_ensembles>0)));
        else
            ECV=1;
            frame_ensembles=1;
            NensemblesOK(CAGindex)=1;
        end
        SAVE_FRAMES{CAGindex}=frame_ensembles;
        ErrorClass(CAGindex)=mean(ECV);
        % MIstFrames(CAGindex)=round(size(Rclust,2)*mean(ECV));
    end
    DelayTime=toc;
    ErrorEnseRation=ErrorClass./NensemblesOK;
    % If Zero->Overffitting or Few Data
    ErrorEnseRation(ErrorEnseRation==0)=1;
    % When 2 ensembles->Missclustering Probably
    ErrorEnseRation(NensemblesOK==2)=1;
    if numel(ErrorEnseRation)>2
        [~,MININIS]=findpeaks(-ErrorEnseRation);
    else
        MININIS=[];
    end
    if isempty(MININIS)
        [~,minErrIndx]=min(ErrorEnseRation);
    else
        minErrIndx=MININIS(1);
    end
    disp('Retrieving Cluster Analysis');
    Rclust=Ractive(:,CAG>=CAGwithAN(minErrIndx));
    % [ActiveCells,~]=find(sum(Rclust,2));
    % Get the Clustered Frames
    % [frame_ensembles]=cluster_analyze(Rclust,SimMethod);
    frame_ensembles=SAVE_FRAMES{minErrIndx}; % saving time!!!!
    signif_frames=OrgigiFrames{minErrIndx}; % saving time!!!!
    % Retrieve Model & ECV -> delete some stuff later...
    NensOK=numel(unique(frame_ensembles(frame_ensembles>0)));
            
    % RE-LABEL Ensembles-> it' got at Ensemble_Sorting Step
    % CAG=sum(R);
    % HebbSequence=Ensembles_Transitions(1,frame_ensembles,signif_frames,CAG,[],0);
    % relabel_frame_ensembles=relabel_ensembles(frame_ensembles(frame_ensembles>0),HebbSequence,'2-freq');
    % frame_ensembles(frame_ensembles>0)=relabel_frame_ensembles;
    [Model,ECV]=Nbayes_Ensembles(R(:,signif_frames(frame_ensembles>0)),frame_ensembles);
    NensOK=numel(unique(frame_ensembles));
    
    % NensOK;
    fprintf('>> Analysis lasted %3.1f seconds  \n',DelayTime);
    fprintf('>> Clustering with %i Ensembles & for %i Coactive Neurons\n',NensOK,CAGwithAN(minErrIndx))
    
    %% Retrieve Ensemble Cluster Array *********************************
    % OUTPUT
    R_Analysis.Data.Data=R';  % Frames x Cells (compatibility wit JP's NN)
    R_Analysis.Peaks.Threshold=CAGwithAN(minErrIndx);
    R_Analysis.Clustering.TotalStates=NensOK;
%     R_Analysis.Clustering.Tree=Tree;
%     R_Analysis.Clustering.Linkage=LinkageMethod;
    R_Analysis.Peaks.Index(CAG>=CAGwithAN(minErrIndx))=1;
    R_Analysis.Clustering.VectorStateIndex=frame_ensembles;  % Sequence of Ensembles: frame_ensembles
    R_Analysis.Classifier.Model=Model;             % Naive Bayes Classifier
    R_Analysis.Classifier.ValidationError=ECV;
else
    % OUTPUT
    R_Analysis.Data.Data=R';  % Frames x Cells (compatibility wit JP's NN)
    R_Analysis.Peaks.Threshold=0;
    R_Analysis.Clustering.TotalStates=0;
%     R_Analysis.Clustering.Tree=Tree;
%     R_Analysis.Clustering.Linkage=LinkageMethod;
    R_Analysis.Peaks.Index=[];
    R_Analysis.Clustering.VectorStateIndex=[];  % Sequence of Ensembles: frame_ensembles
    R_Analysis.Classifier.Model=[];             % Naive Bayes Classifier
    R_Analysis.Classifier.ValidationError=[];
    CAGwithAN=[];
end
% fprintf('>> Clustering with %i Ensembles \n& for %i Coactive Neurons\n',Nensembles,CAG_TH)
% fprintf('\n>>Algorithm has found Neural Ensembles.\n')
%% FAST PLOTTING *******************************************************
fprintf('\n>>Plotting Ensembles:')
ImageEnsembles(R_Analysis,1);
if ~isempty(CAGwithAN)
    signif_frames=find(CAG>=CAGwithAN(minErrIndx));
    [New_Order_Clustering,~]=OrderClusters(frame_ensembles(frame_ensembles>0),...
        signif_frames(frame_ensembles>0),R',NensOK);
    R_Sorted_Analysis=R_Analysis;
    R_Sorted_Analysis.Data.Data=R(New_Order_Clustering,:)';
    ImageEnsembles(R_Sorted_Analysis,0);
    % Set Indexes at sorted Ensembles
    figHandles = findobj('Type', 'figure');
    figHandles(3).Children(2).YTick=1:max(New_Order_Clustering);
    figHandles(3).Children(2).YTickLabel=New_Order_Clustering;
    fprintf(' ready!\n')
else
    
    fprintf('empty\n')
end

fprintf('\n>>Ensembles Done\n')