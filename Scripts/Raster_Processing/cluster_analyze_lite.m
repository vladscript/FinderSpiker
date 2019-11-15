% Function that analyze if the two main branches of the
% cluster analysis (using hierarchichal clustering)
% have other cluster that helps to classfify better
% the neural ensembles
% Input
%   Rclust:     Thresholded Raster to Analyze
%   SimMethod:  Simmilarity Method: 'hamming'
% Ouput
% frame_ensembles Labeled Frames of the Neural Ensemble
function [frame_ensembles]=cluster_analyze_lite(Rclust)
%% Setup
Load_Default_Clustering;
% BINARY CLASSIFICATION ***************************************
% Algorithm for computing distance between clusters.
% see: >>doc linkage
%% Initial Tree ************************************************
Distance=squareform(pdist(Rclust',SimMethod));
Sim=1-Distance;
LinkageMethod=HBestTree_JPplus(Sim);    % Output
Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
%% Setup Loop
Nensambles=2;
MoreEnsembles=true;
preDOMI=[];
preActiveFrames=[];
[ActiveCellsTotal,~]=find(sum(Rclust,2));
StopbyOverthrown=false;
StopbyRedundant=false;
% Main Loop *****************************************************
while MoreEnsembles
    AllCellEns=[];
    NeuroVectors=zeros(numel(ActiveCellsTotal),Nensambles);
    EnsembledNeurons={};
    fprintf('>>Trying %i Neural Ensembles\n',Nensambles)
    % Hierarchichal Clustering:
    frame_ensembles=cluster(Tree,'maxclust',Nensambles); % START
    %  Checking Separation of the Most Frequent:
    tbl=tabulate(frame_ensembles);
    [ActiveFrames,DOMI]=max(tbl(:,2));
    if Nensambles>2
        % ask if it is the same ENS DOMINANT
        if ~ismember(DOMI,preDOMI);
            MoreEnsembles=false;
            StopbyOverthrown=true;
            fprintf('>>Overthrown Ensemble.\n')
        end
    end
    preDOMI=[preDOMI;DOMI];
    preActiveFrames=[preActiveFrames;ActiveFrames];
    disp([preDOMI';preActiveFrames'])
    
    if MoreEnsembles
        % Checking Redundant Ensembles
        for nn=1:Nensambles
            EnsembledNeurons{nn}=find(sum(Rclust(:,frame_ensembles==nn),2));
            NeuroVectors(EnsembledNeurons{nn},nn)=1;
        end
        EnsemblesIndx=isinensemble(NeuroVectors);
        % All-Cells Ensemble:
        AllCellEns=find(sum(NeuroVectors)==size(Rclust,1));
        if ~isempty(EnsemblesIndx)
            disp('>>Redundant Ensemble(s) Found!!!')
            [NRepEns,~]=size(EnsemblesIndx);
            FracEns=zeros(NRepEns,1);
            for nE=1:NRepEns
                fprintf('\nEnsemble: %i is in Ensemble %i ',EnsemblesIndx(nE,1),EnsemblesIndx(nE,2));
                FractionofEnsemble=numel(EnsembledNeurons{EnsemblesIndx(nE,1)})/numel(EnsembledNeurons{EnsemblesIndx(nE,2)});
                FracEns(nE)=FractionofEnsemble>parseq;
                fprintf('a %3.2f %% \n',100*FractionofEnsemble);
            end
            if isempty(AllCellEns) && sum(FracEns)>0
                % Stop Only if there is NON all-cell Ensemble
                MoreEnsembles=false;
                StopbyRedundant=true;
                fprintf('>>No more Ensmebles.')
            elseif ~isempty(AllCellEns)
                fprintf('>>Ensemble %i contains all cells\n',AllCellEns)
            else
                fprintf('>>Redundant but Separable\n')
            end
        end
    end
    Nensambles=Nensambles+1;
end
% Setting Number of Ensembles:

if StopbyOverthrown
    NensamblesOK=Nensambles-1;
else
    NensamblesOK=Nensambles-2;
end
if NensamblesOK<2
    NensamblesOK=2;
end
fprintf('>>Set to %i Neural Ensembles',NensamblesOK)
frame_ensembles=cluster(Tree,'maxclust',NensamblesOK); % START
fprintf(' Ready.\n')


% 
% 