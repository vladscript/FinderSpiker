% Function that analyze if the two main branches of the
% cluster analysis (using hierarchichal clustering)
% have other cluster that helps to classfify better
% the neural ensembles
% Input
%   Rclust: Raster to analyze
%   SimMethod: Simmilarity Method
% Ouput
% frame_ensembles Labeled Frames of the Neural Ensemble
function [frame_ensembles]=cluster_analyze(Rclust,SimMethod)
    % BINARY CLASSIFICATION ***************************************
    Distance=squareform(pdist(Rclust',SimMethod));
    Sim=1-Distance;
    LinkageMethod=HBestTree_JPplus(Sim);    % Output
    Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
    frame_ensembles=cluster(Tree,'maxclust',2); % START
    % Analyze the Main two Tree Branches  #################################
    TotalEns=0; SUB_ENS_frames={}; DhammEns={};
    for n=1:2
        % Intra Ensembles
        Rclustens=Rclust(TwoEnsembledNeurons{n},frame_ensembles==1+TotalEns);
        [ActiveCellsEns,~]=find(sum(Rclustens,2));
        DistanceEns=squareform(pdist(Rclustens',SimMethod));
        SimEns=1-DistanceEns;
        LinkageMethodEns=HBestTree_JPplus(SimEns);    % Output
        TreeEns=linkage(squareform(DistanceEns,'tovector'),LinkageMethodEns);
        % Initialize stuff
        EnsambleOK=true; Nensembles=2; preECVens=1;
        while EnsambleOK
            fprintf('>>> Clustering for  %i Ensembles \n',Nensembles);
            NeuroVectors=zeros(numel(ActiveCellsEns),Nensembles);
            EnsembledNeurons={};
            SUB_frame_ensembles=cluster(TreeEns,'maxclust',Nensembles); % START
            [~,ECVens]=Nbayes_Ensembles(Rclustens,SUB_frame_ensembles);

            for nn=1:Nensembles
                EnsembledNeurons{nn}=find(sum(Rclustens(:,SUB_frame_ensembles==nn),2));
                NeuroVectors(EnsembledNeurons{nn},nn)=1;
            end

            % Redundant Ensembles
            EnsemblesIndx=isinensemble(NeuroVectors);
            if ~isempty(EnsemblesIndx)
                disp('>>Redundant Ensemble(s) Found!!!')
                [NRepEns,~]=size(EnsemblesIndx);
                for nE=1:NRepEns
                    fprintf('\nEnsemble: %i is in Ensemble %i\n',EnsemblesIndx(nE,1),EnsemblesIndx(nE,2));
                end
                EnsambleOK=false;
            end
            disp('>>Done.')

            % 1-neuron Ensembles
            if sum(sum(NeuroVectors)==1)>0
                disp('1-neuron Ensembles');
                EnsambleOK=false;
            end

            % Increased Error Classification
            if preECVens-ECVens<0
                EnsambleOK=false;
            end
            preECVens=ECVens;
            if EnsambleOK
                % Only increas if it's OK
                Nensembles=Nensembles+1;
            end
        end
        Nensembles=Nensembles-1;
        if Nensembles>1
            SUB_ENS_frames{n}=cluster(TreeEns,'maxclust',Nensembles)+TotalEns; % START
            EnsembledNeurons={};
            NeuroVectors=zeros(numel(ActiveCellsEns),Nensembles);
            for nn=1:Nensembles
                EnsembledNeurons{nn}=find(sum(Rclustens(:,SUB_ENS_frames{n}==nn+TotalEns),2));
                NeuroVectors(EnsembledNeurons{nn},nn)=1;
            end
            % INTRA CLUSTER DISTANCE
            DhammEns{n}=pdist(NeuroVectors',SimMethod);
        else
            disp('<<No sub-ensambles found>>')
            SUB_ENS_frames{n}=frame_ensembles(frame_ensembles==TotalEns+1);
        end
        frame_ensembles(frame_ensembles>n+TotalEns)=frame_ensembles(frame_ensembles>n+TotalEns)+Nensembles-1; 
        frame_ensembles(frame_ensembles==TotalEns+1)=SUB_ENS_frames{n};
        TotalEns=TotalEns+Nensembles;

    end