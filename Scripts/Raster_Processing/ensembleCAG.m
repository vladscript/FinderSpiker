% Function that ignores low-active cells in Ensembles
% Input
%   frame_ensemble: sequence of ensembles
%   Rclust:         Matriz Activity above CAG threshold
% Output
%   reframe_ensembles: zeros where low-active cells fired
function reframe_ensembles=ensembleCAG(frame_ensemble,Rclust)
reframe_ensembles=frame_ensemble;
disp('>>Anlayzing CAG ensemble')
EnsmblesList=unique(frame_ensemble);
Nensembles=numel(EnsmblesList);
SimMethod='hamming';
for n=1:Nensembles
    FramesEnsemble=find(frame_ensemble==EnsmblesList(n));
    Rcag=Rclust(:,frame_ensemble==EnsmblesList(n));
    EnsembleVector=zeros(size(Rcag,1),1);
    EnsembleVector(sum(Rcag,2)>0)=1;
    Rensembled=Rcag(sum(Rcag,2)>0,:);
    RoAensembled=sum(Rensembled,2)/size(Rensembled,2);
    % Get Most ACtive Neurons During Ensemble Instances
    % we expect high RoA due to it's only Enemble Activity
    Dlist=pdist(Rensembled',SimMethod);
    preCut=max(Dlist);
    % Apply Threshold RoA
    if median(RoAensembled)<0.25
        ThRoA=median(RoAensembled);
    else
        ThRoA=0.25;
    end
    
    RejCEns=find(EnsembleVector);
    RejectedCells=RejCEns(RoAensembled<ThRoA);
    %RejectFrames=find(sum(Rcag(RejectedCells,:)));
    reframe_ensembles(FramesEnsemble(sum(Rcag(RejectedCells,:))>0))=0;
%     % Checking:
%     [~,SortIndex]=sort(RoAensembled)
%     Plot_Raster_Ensembles(Rensembled(SortIndex,:))
%     disp(ThRoA)
%     RejectCells=find(RoAensembled<ThRoA);
%     disp(RejectCells)
%     pause;
end
% Check Final result
% Plot_Raster_Ensembles(Rclust(:,reframe_ensembles>0))
