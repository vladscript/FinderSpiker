


%% pre - VERSION

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