% Ensembles Transitions: Hebbian Pathway
% Input
%   signif_frames:  Frames with Significative CoACtivity
%   labels_frames:  Ensemble Lalbel
%   ColorState:     Ensemble Color
%   If plot=1 else 0
% Output
%   ensemble_index: Sequence of Ensembles                
function [ensemble_index]=Ensembles_Transitions(fs,labels_frames,signif_frames,ColorState,ifplot)
frames_ensembles=[];
aux=1;
time_ensemble=[];
ensemble_index=[];
labels_frames=[labels_frames;1000];

for i=1:length(signif_frames);
    if labels_frames(i+1)==labels_frames(i)
        frames_ensembles=[frames_ensembles,signif_frames(i)];
    else
        if isempty(frames_ensembles)
            frames_ensembles=signif_frames(i);
        end
        time_ensemble(aux,1)=round(mean(frames_ensembles));
        ensemble_index(aux,1)=labels_frames(i);
        aux=aux+1;
        frames_ensembles=[];
    end 
end

if ifplot
    %% Plot in Minutes
    Axis_details=gca;
    fig_Ensembles_Transitions=figure('Position', [415 172 560 122],...
        'Name','Hebbian Pathways');
    plot(time_ensemble/fs/60,ensemble_index,'k','LineWidth',3); hold on
    Ntran=length(ensemble_index);
    for i=1:Ntran
        plot(time_ensemble(i)/fs/60,ensemble_index(i),'o',...
            'MarkerEdgeColor',ColorState(ensemble_index(i),:),...
            'MarkerFaceColor',ColorState(ensemble_index(i),:),...
            'MarkerSize',13); 
    end
    hold off;
    axis([Axis_details.XLim,1,max(ensemble_index)])
end