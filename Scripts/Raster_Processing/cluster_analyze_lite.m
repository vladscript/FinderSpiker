% Function that analyze if the two main branches of the
% cluster analysis (using hierarchichal clustering)
% have other cluster that helps to classfify better
% the neural ensembles
% Input
%   Rclust:     Thresholded Raster to Analyze
%   SimMethod:  Simmilarity Method: 'hamming'
% Ouput
% frame_ensembles Labeled Frames of the Neural Ensemble
function [frame_ensembles]=cluster_analyze_lite(Rclust,SimMethod)
    % BINARY CLASSIFICATION ***************************************
    Distance=squareform(pdist(Rclust',SimMethod));
    Sim=1-Distance;
    LinkageMethod=HBestTree_JPplus(Sim);    % Output
    Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
    frame_ensembles=cluster(Tree,'maxclust',2); % START
    % Check if there single frame ensembles and ignore them
    tbl=tabulate(frame_ensembles);
    SingleEns=find(tbl(:,2)==1);
    Rnewclust=Rclust;
    aux=1; % frames2ignore=[];
    while ~isempty(SingleEns)
        fprintf('>>Ignoring Single Frame Ensemble %i-th trial\n',aux)
        NsingEns=numel(SingleEns);
        frames2ignore=[];
        for k=1:NsingEns
            frames2ignore=[frames2ignore,find(frame_ensembles==SingleEns(k))];
        end
        Rnewclust(:,frames2ignore)=0; % Ignore these vector
        Distance=squareform(pdist(Rnewclust',SimMethod));
        Sim=1-Distance;
        LinkageMethod=HBestTree_JPplus(Sim);    % Output
        Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
        frame_ensembles=cluster(Tree,'maxclust',2); % START
        % Check if there single frame ensembles and ignore them
        tbl=tabulate(frame_ensembles);
        SingleEns=find(tbl(:,2)==1);
        aux=aux+1;
    end
    % frame_ensembles(frames2ignore)=0;
    
    % NeuroVectors=zeros(numel(ActiveCells),Nensembles);
    TwoEnsembledNeurons={};
    for nn=1:2
        TwoEnsembledNeurons{nn}=find(sum(Rclust(:,frame_ensembles==nn),2));
        %TwoNeuroVectors(TwoEnsembledNeurons{nn},nn)=1;
    end
    % Analyze the Main two Tree Branches  #################################
    TotalEns=0; SUB_ENS_frames={}; DhammEns={};
    for n=1:2
        % Clusterize only the biggest ensemble ONLY
        if numel(frame_ensembles(frame_ensembles==1+TotalEns))<numel(frame_ensembles)/2
            EnsambleOK=false; Nensembles=2; preMAXens=numel(frame_ensembles);
        else
            % Intra Ensembles
            Rclustens=Rclust(TwoEnsembledNeurons{n},frame_ensembles==1+TotalEns);
            % Identify Frames Where Vector are not Very Simmilar
            % [~,MaxDistanceIV]=nonsimilarframes(Rclustens,SimMethod,0.5);
            % NeuralVectorNonSim=find(sum(Rclustens(:,NonSimilarIV),2))

            [ActiveCellsEns,~]=find(sum(Rclustens,2));
            DistanceEns=squareform(pdist(Rclustens',SimMethod));
            SimEns=1-DistanceEns;
            LinkageMethodEns=HBestTree_JPplus(SimEns);    % Output
            TreeEns=linkage(squareform(DistanceEns,'tovector'),LinkageMethodEns);
            % Initialize stuff
            EnsambleOK=true; Nensembles=2; preMAXens=numel(frame_ensembles);
        end
        while EnsambleOK
            fprintf('>>> Clustering for  %i Ensembles \n',Nensembles);
            NeuroVectors=zeros(numel(ActiveCellsEns),Nensembles);
            EnsembledNeurons={};
            SUB_frame_ensembles=cluster(TreeEns,'maxclust',Nensembles); % START
            
            tbl=tabulate(SUB_frame_ensembles);
            SingleEns=find(tbl(:,2)<2);
            Rnewclust=Rclustens;
            aux=1; % frames2ignore=[];
            while and(~isempty(SingleEns),aux<=5)
                fprintf('>>Ignoring Single Frame Ensemble %i-th trial\n',aux)
                NsingEns=numel(SingleEns);
                frames2ignore=[];
                for k=1:NsingEns
                    frames2ignore=[frames2ignore,find(SUB_frame_ensembles==SingleEns(k))];
                end
                Rnewclust(:,frames2ignore)=0; % Ignore these vectors
                Distance=squareform(pdist(Rnewclust',SimMethod));
                Sim=1-Distance;
                LinkageMethod=HBestTree_JPplus(Sim);    % Output
                Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
                SUB_frame_ensembles=cluster(Tree,'maxclust',Nensembles); % START
                % Check if there single frame ensembles and ignore them
                tbl=tabulate(SUB_frame_ensembles);
                SingleEns=find(tbl(:,2)<2);
                aux=aux+1;
            end
            
            % Single Frame Ensembles
            if ~isempty(SingleEns)
                EnsambleOK=false;
                disp('>> Unique-frame Ensemble: rejected.')
            end
            
            for nn=1:Nensembles
                EnsembledNeurons{nn}=find(sum(Rclustens(:,SUB_frame_ensembles==nn),2));
                NeuroVectors(EnsembledNeurons{nn},nn)=1;
            end
            
            
            % 1-neuron Ensembles
            if sum(sum(NeuroVectors)==1)>0
                disp('1-neuron Ensembles');
                EnsambleOK=false;
            end

            % Decreas Max Counts Ensemble
            tbl=tabulate(SUB_frame_ensembles);
            MaxEns= max(tbl(:,2));
            if preMAXens-MaxEns<0 || max(tbl(:,3))<50
                EnsambleOK=false;
            end
            preMAXens=MaxEns;
            
            if EnsambleOK
                % It only increases if it's OK
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
    %% Work-arounds noIA:ifs and more ifs**********************************
    tbl=tabulate(frame_ensembles); % Get counts for each ensemble
    % Single-frame Ensmebles ##############################################
    % may cause troubla t saving time features of ensembles
    SingEns=find(tbl(:,2)<2);
    if isempty(SingEns)
        disp('>>No Single-Frame Ensembles left')
    else
        disp('>>Adjusting Single-Frame Ensembles left')
        frame_ensembles_copy=frame_ensembles;
        EnsembledNeurons={};
        NeuroVectors=zeros(size(Rclust,1),TotalEns);
        for nn=1:TotalEns
            EnsembledNeurons{nn}=find(sum(Rclust(:,frame_ensembles==nn),2));
            NeuroVectors(EnsembledNeurons{nn},nn)=1;
        end
        DistEns= squareform( pdist(NeuroVectors',SimMethod) );
        Nsing=numel(SingEns);
        for n=1:Nsing
            [~,EqEns]=max(DistEns(SingEns(n),:));
            frame_ensembles_copy(frame_ensembles==SingEns(n))=EqEns;
            fprintf('>>Ensemble %i-> Assigned to ->Ensemble %i \n',SingEns(n),EqEns)
        end
        % Re-adjust frame labels
        labens=unique(frame_ensembles_copy);
        for nn=1:numel(labens)
            frame_ensembles(frame_ensembles_copy==labens(nn))=nn;
        end
    end
    %% END
    disp('>>Analysis Completed.')