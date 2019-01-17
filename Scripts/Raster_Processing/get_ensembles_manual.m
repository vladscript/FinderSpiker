% Function that make Hierarchichal Cluster 
% from Matrix Activity
% Input
%   Raster: Matrix Cells x Frames
%   CAG Threshold: Scalar
%   Simmilaarity Method : 'hamming', etc
% Ouput
% R_Analysis Output Strcuture Output
function R_Analysis = get_ensembles_manual(R,CAGth,SimMethod)
    ActiveCells= sum(R,2);
    % GET ONLY ACTIVE CELLS!
    Ractive=R(ActiveCells>0,:)';
    CAG=sum(Ractive);
    
    MaxCAG=max(CAG);

    if MaxCAG<CAGth
        fprintf('\n>> Raster with LOW Co-Activity \n\n');
        R_Analysis.Data.Data=R'; % Frames x Cells (compatibility wit JP's NN)
        fprintf('>> Finishing script...');
        Analyze=false;
    else
        Analyze=true;
    end
    if Analyze
        disp('>>Analysing...')
        Rclust=Ractive(:,CAG>=CAGth);
        % Get the Clustered Frames
        [frame_ensembles]=cluster_analyze(Rclust,SimMethod);
        NensemblesOK=numel(unique(frame_ensembles));
        % RE-LABEL Ensembles
        AppearSequence=unique(frame_ensembles,'stable');
        relabel_frame_ensembles=zeros(size(frame_ensembles));
        for n=1:NensemblesOK
            relabel_frame_ensembles(frame_ensembles==AppearSequence(n))=n;
        end
        [Model,ECV]=Nbayes_Ensembles(Rclust,relabel_frame_ensembles);
        
        % Save Stuff
        R_Analysis.Data.Data=R;  % Frames x Cells (compatibility wit JP's NN)
        R_Analysis.Peaks.Threshold=CAGth;
        R_Analysis.Clustering.TotalStates=NensemblesOK;
        % R_Analysis.Clustering.Tree=Tree;
        % R_Analysis.Clustering.Linkage=LinkageMethod;
        R_Analysis.Peaks.Index(CAG>=CAGth)=1;
        R_Analysis.Clustering.VectorStateIndex=relabel_frame_ensembles;
        R_Analysis.Classifier.Model=Model;
        R_Analysis.Classifier.ValidationError=ECV;
        % PLOT ENSEMBLES UNSORT & RESORTED
        fprintf('\n>>Plotting Ensembles:')
        ImageEnsembles(R_Analysis);
        signif_frames=find(CAG>=CAGth);
        [New_Order_Clustering,~]=OrderClusters(relabel_frame_ensembles,signif_frames,R,NensemblesOK);
        R_Sorted_Analysis=R_Analysis;
        R_Sorted_Analysis.Data.Data=R(:,New_Order_Clustering);
        ImageEnsembles(R_Sorted_Analysis);
        fprintf('\n>>Done')
    end
    disp('>>Done.')
end