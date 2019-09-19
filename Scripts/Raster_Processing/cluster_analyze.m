% Function that analyze if the two main branches of the
% cluster analysis (using hierarchichal clustering)
% have other cluster that helps to classfify better
% the neural ensembles
% Input
%   Rclust:         Thresholded Raster to Analyze
% Ouput
%   frame_ensembles   Labeled Frames of the Neural Ensemble
function [frame_ensembles_ok]=cluster_analyze(Rclust)
%% Setup ********************************************
Load_Default_Clustering;
%% First Clustering: ********************************
[frame_ensembles]=cluster_analyze_lite(Rclust);

%% RE-FRAME *****************************************
% reframe_ensembles=ensembleCAG(frame_ensembles,Rclust);
% disp('ass')
FramesOut=find(frame_ensembles>0); % From Input Raster--------><
% Rnewclust=Rclust(:,reframe_ensembles>0);
% ActNeu=find(sum(Rnewclust,2));
% re set raster and frames:
% Rnewclust=Rnewclust(ActNeu,:);
% frame_ensembles=reframe_ensembles(reframe_ensembles>0);
CAG=sum(Rclust);

%% Joint Simmilar Ensembles
fprintf('>>Checking Low-Frequent Ensemble')
% Make it Recurrent as ensmebles assembly together

% To enter the Loop
hebbtbl=tabulate(frame_ensembles);
frame_ensembles_copy=frame_ensembles;
smallens=find(hebbtbl(:,3)<RatioHebb);
% Neural Vectors and Distances Intra Ensembles
NensOK=numel(unique(frame_ensembles));
NeuroVectors=zeros(size(Rclust,1),NensOK);
for nn=1:NensOK
    NeuroVectors(sum(Rclust(:,frame_ensembles==nn),2)>0,nn)=1;
end
DistIV=squareform(pdist(NeuroVectors',SimMethod));
SimIV=1-DistIV;
while ~isempty(smallens)
    ActualEns=smallens(1);
    OtherEnsembles=setdiff(1:NensOK,ActualEns);
    [~,EnseRexIndx]=max(SimIV(ActualEns,OtherEnsembles));
    EnsRx=OtherEnsembles(EnseRexIndx);
    % Join:
    fprintf('>>Joining: Ensemble %i & %i\n',ActualEns,EnsRx);
    MergeingEnsembles=[min(ActualEns,EnsRx),max(ActualEns,EnsRx)];
    % Update Ensembles
    frame_ensembles(frame_ensembles_copy==MergeingEnsembles(2))=MergeingEnsembles(1);
    frame_ensembles_copy=frame_ensembles;
    % Avoid Enemble with Zero vectors:
    if MergeingEnsembles(2)<NensOK
        Nshifts=NensOK-MergeingEnsembles(2);
        for nn=1:Nshifts
            frame_ensembles(frame_ensembles_copy==MergeingEnsembles(2)+nn)=MergeingEnsembles(2)+nn-1;
        end
    end
    frame_ensembles_copy=frame_ensembles;
    % Re -get table:
    hebbtbl=tabulate(frame_ensembles);
    smallens=find(hebbtbl(:,3)<RatioHebb);
    % Neural Vectors and Distances Intra Ensembles
    NensOK=numel(unique(frame_ensembles));
    NeuroVectors=zeros(size(Rclust,1),NensOK);
    for nn=1:NensOK
        NeuroVectors(sum(Rclust(:,frame_ensembles==nn),2)>0,nn)=1;
    end
    DistIV=squareform(pdist(NeuroVectors',SimMethod));
    SimIV=1-DistIV;
end
fprintf(': Done.\n')

%% WORK AROUND ********************************************************
% Excesive Dominant Ensemble ##########################################
% This causes *cluster artifacts*
fprintf('>>Checking Super Ensemble')
tbl=tabulate(frame_ensembles); % Get counts for each ensemble
Nenspre=numel(unique(frame_ensembles));
framesCAG=find(CAG>=min(CAG));
[TimePer,ExcEns]=max(tbl(:,3)); % Percentage of presence
if TimePer>TimeDommEns
    disp('>>Possible Extremely Dominant Ensemble')
    % Is it alterned or sequentlty activated
    samples=framesCAG(frame_ensembles==ExcEns);
    diffsamples=diff(samples);
    if numel(diffsamples(diffsamples==1))>numel(samples)*parseq
        fprintf('>>Ensemble Dominant Permanently Activated %3.1f %% of the time\n',100*numel(diffsamples(diffsamples==1))/numel(samples))
        % Re-cluster Dominant Ensemble
        RclustDom=Rclust(:,samples);
        CAGclusDom=sum(RclustDom);
        [re_frame_ensembles]=cluster_analyze_lite(RclustDom);
        hebbtbl=tabulate(re_frame_ensembles);
        smallens=find(hebbtbl(:,3)<RatioHebb);
        for k=1:numel(smallens)
            re_frame_ensembles(re_frame_ensembles==smallens(k))=smallens(1);
        end
        Ensd=unique(re_frame_ensembles);
        kens=numel(Ensd);
        re_frame_ensembles_copy=re_frame_ensembles;
        for k=1:kens
            re_frame_ensembles(re_frame_ensembles_copy==Ensd(k))=k;
        end
        Nensoffset=numel(unique(re_frame_ensembles));
        frame_ensembles_copy=frame_ensembles+Nensoffset;
        frame_ensembles_copy(frame_ensembles_copy==ExcEns+Nensoffset)=re_frame_ensembles;
        frame_ensembles=frame_ensembles_copy;
        Ensd=unique(frame_ensembles);
        kens=numel(Ensd);
        frame_ensembles_copy=frame_ensembles;
        for k=1:kens
            frame_ensembles(frame_ensembles_copy==Ensd(k))=k;
        end
    else
        disp('>>Ensemble Dominant but Alternating.')

    end
else
    disp('>>Non-Dominant Ensemble Found.')
end
fprintf('>>Done\n\n');
%% OUTPUT 
frame_ensembles_ok=zeros(size(Rclust,2),1);
frame_ensembles_ok(FramesOut)=frame_ensembles;

% Check if there single frame ensembles and ignore them

% SingleEns=find(tbl(:,2)==1);    % Outlier Ensemble(s)
% Rnewclust=Rclust;
% aux=1; % frames2ignore=[];
% while and(~isempty(SingleEns),aux<MaxEnsembles)
%     fprintf('>>Ignoring Single Frame Ensemble %i-th trial\n',aux)
%     NsingEns=numel(SingleEns);
%     frames2ignore=[];
%     for k=1:NsingEns
%         frames2ignore=[frames2ignore,find(frame_ensembles==SingleEns(k))];
%     end
%     Rnewclust(:,frames2ignore)=0; % Ignore these vector
%     Distance=squareform(pdist(Rnewclust',SimMethod));
%     Sim=1-Distance;
%     LinkageMethod=HBestTree_JPplus(Sim);    % Output
%     Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
%     frame_ensembles=cluster(Tree,'maxclust',2); % START
%     % Check if there single frame ensembles and ignore them
%     tbl=tabulate(frame_ensembles);
%     SingleEns=find(tbl(:,2)==1);
%     aux=aux+1;
% end
% if ~isempty(SingleEns)
%     disp('>>WARNING: single-frame Ensembles')
% end
% % frame_ensembles(frames2ignore)=0;
% 
% % NeuroVectors=zeros(numel(ActiveCells),Nensembles);
% TwoEnsembledNeurons={};
% for nn=1:2
%     TwoEnsembledNeurons{nn}=find(sum(Rclust(:,frame_ensembles==nn),2));
%     CellRatioEns(nn)=numel(TwoEnsembledNeurons{nn})/size(Rclust,1);
% end
% % CoreCell=intersect(TwoEnsembledNeurons{1},TwoEnsembledNeurons{2});
% % Analyze the Main two Tree Branches  #################################
% TotalEns=0; SUB_ENS_frames={}; DhammEns={};
% for n=1:2
%     % Intra Ensembles
%     Rclustens=Rclust(TwoEnsembledNeurons{n},frame_ensembles==1+TotalEns);
%     if size(Rclustens,2)==1
%         EnsambleOK=false;
%         Nensembles=2;
%     else
%         [ActiveCellsEns,~]=find(sum(Rclustens,2));
%         DistanceEns=squareform(pdist(Rclustens',SimMethod));
%         SimEns=1-DistanceEns;
%         LinkageMethodEns=HBestTree_JPplus(SimEns);    % Output
%         TreeEns=linkage(squareform(DistanceEns,'tovector'),LinkageMethodEns);
%         % Initialize stuff
%         EnsambleOK=true; Nensembles=2; preECVens=1;
%     end
%     % Identify Frames Where Vector are not Very Simmilar
%     % [~,MaxDistanceIV]=nonsimilarframes(Rclustens,SimMethod,0.5);
%     % NeuralVectorNonSim=find(sum(Rclustens(:,NonSimilarIV),2))
% 
%     while EnsambleOK
%         fprintf('>>> Clustering for  %i Ensembles \n',Nensembles);
%         NeuroVectors=zeros(numel(ActiveCellsEns),Nensembles);
%         EnsembledNeurons={};
%         SUB_frame_ensembles=cluster(TreeEns,'maxclust',Nensembles); % START
% 
%         tbl=tabulate(SUB_frame_ensembles);
%         SingleEns=find(tbl(:,2)<2);
%         Rnewclust=Rclustens;
%         aux=1; % frames2ignore=[];
%         while and(~isempty(SingleEns),aux<=5)
%             fprintf('>>Ignoring Single Frame Ensemble %i-th trial\n',aux)
%             NsingEns=numel(SingleEns);
%             frames2ignore=[];
%             for k=1:NsingEns
%                 frames2ignore=[frames2ignore,find(SUB_frame_ensembles==SingleEns(k))];
%             end
%             Rnewclust(:,frames2ignore)=0; % Ignore these vectors
%             Distance=squareform(pdist(Rnewclust',SimMethod));
%             Sim=1-Distance;
%             LinkageMethod=HBestTree_JPplus(Sim);    % Output
%             Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
%             SUB_frame_ensembles=cluster(Tree,'maxclust',Nensembles); % START
%             % Check if there single frame ensembles and ignore them
%             tbl=tabulate(SUB_frame_ensembles);
%             SingleEns=find(tbl(:,2)<2);
%             aux=aux+1;
%         end
% 
% 
%         [~,ECVens]=Nbayes_Ensembles(Rclustens,SUB_frame_ensembles);
%         for nn=1:Nensembles
%             EnsembledNeurons{nn}=find(sum(Rclustens(:,SUB_frame_ensembles==nn),2));
%             NeuroVectors(EnsembledNeurons{nn},nn)=1;
%         end
%         % Single Frame Ensembles
%         if ~isempty(SingleEns)
%             EnsambleOK=false;
%             disp('>> Unique-frame Ensemble: rejected.')
%         end
% 
%         % Redundant Ensembles
%         EnsemblesIndx=isinensemble(NeuroVectors);
%         % All-Cells Ensemble:
%         AllCellEns=find(sum(NeuroVectors)==size(Rclustens,1));
%         if ~isempty(EnsemblesIndx)
%             disp('>>Redundant Ensemble(s) Found!!!')
%             [NRepEns,~]=size(EnsemblesIndx);
%             for nE=1:NRepEns
%                 fprintf('\nEnsemble: %i is in Ensemble %i\n',EnsemblesIndx(nE,1),EnsemblesIndx(nE,2));
%             end
%             if isempty(AllCellEns)
%                 % Stop Only if there is no all-cell ensembles
%                 EnsambleOK=false;
%             end
%         end
%         disp('>>Done.')
% 
%         % 1-neuron Ensembles (no way to happen, CAG>1)
%         if sum(sum(NeuroVectors)==1)>0
%             disp('1-neuron Ensembles');
%             EnsambleOK=false;
%         end
% 
%         % Increased Error Classification
%         if preECVens-ECVens<0
%             disp('>> Worsening Classification');
%             EnsambleOK=false;
%         end
%         preECVens=ECVens;
%         if EnsambleOK
%             % Only increas if it's OK
%             Nensembles=Nensembles+1;
%         end
%     end
%     Nensembles=Nensembles-1;
%     if Nensembles>1
%         SUB_ENS_frames{n}=cluster(TreeEns,'maxclust',Nensembles)+TotalEns; % START
%         EnsembledNeurons={};
%         NeuroVectors=zeros(numel(ActiveCellsEns),Nensembles);
%         for nn=1:Nensembles
%             EnsembledNeurons{nn}=find(sum(Rclustens(:,SUB_ENS_frames{n}==nn+TotalEns),2));
%             NeuroVectors(EnsembledNeurons{nn},nn)=1;
%         end
%         % INTRA CLUSTER DISTANCE
%         DhammEns{n}=pdist(NeuroVectors',SimMethod);
%     else
%         disp('<<No sub-ensambles found>>')
%         SUB_ENS_frames{n}=frame_ensembles(frame_ensembles==TotalEns+1);
%     end
%     frame_ensembles(frame_ensembles>n+TotalEns)=frame_ensembles(frame_ensembles>n+TotalEns)+Nensembles-1; 
%     frame_ensembles(frame_ensembles==TotalEns+1)=SUB_ENS_frames{n};
%     TotalEns=TotalEns+Nensembles;
% 
% end
% %% Work-arounds noIA:ifs and more ifs**********************************
% tbl=tabulate(frame_ensembles); % Get counts for each ensemble
% % Single-frame Ensmebles ##############################################
% % may cause troubla t saving time features of ensembles
% SingEns=find(tbl(:,2)<2);
% if isempty(SingEns)
%     disp('>>No Single-Frame Ensembles left')
% else
%     disp('>>Adjusting Single-Frame Ensembles left')
%     frame_ensembles_copy=frame_ensembles;
%     EnsembledNeurons={};
%     NeuroVectors=zeros(size(Rclust,1),TotalEns);
%     for nn=1:TotalEns
%         EnsembledNeurons{nn}=find(sum(Rclust(:,frame_ensembles==nn),2));
%         NeuroVectors(EnsembledNeurons{nn},nn)=1;
%     end
%     DistEns= squareform( pdist(NeuroVectors',SimMethod) );
%     Nsing=numel(SingEns);
%     for n=1:Nsing
%         [~,EqEns]=max(DistEns(SingEns(n),:));
%         frame_ensembles_copy(frame_ensembles==SingEns(n))=EqEns;
%         fprintf('>>Ensemble %i-> Assigned to ->Ensemble %i \n',SingEns(n),EqEns)
%     end
%     % Re-adjust frame labels
%     labens=unique(frame_ensembles_copy);
%     for nn=1:numel(labens)
%         frame_ensembles(frame_ensembles_copy==labens(nn))=nn;
%     end
% end
% %% Anayze CAG Ensemble Cell Variability
% reframe_ensembles=ensembleCAG(frame_ensembles,Rclust);
% frame_ensembles=reframe_ensembles;
% %% END
% disp('>>Analysis Completed.')