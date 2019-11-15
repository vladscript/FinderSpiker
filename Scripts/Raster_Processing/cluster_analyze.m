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
if ismember(0,unique(frame_ensembles))
    disp('>>Cluster Zero<<')
    frame_ensembles=frame_ensembles+1;
    %Nensambles=numel(unique(frame_ensembles))
end
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

% 