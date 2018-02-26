%% Plotting States Colors  *********************************************
% Request
% Plot Raster
% Input
%   Experiment:    Raster Concatenated and Ensemble Sorted
%   labels_frames: Ensemble Labels in each fram
%   signif_frames: Significative Coactivity Frames
%   ColorState:    Color of The Ensembles
%   THR: Coactivity Threshold
% Ouput
% Plot on the Raster plot
function Plot_State_Colors(labels_frames,signif_frames,ColorState,Experiment,THR,fs,CoAc,indexes)
% Search of figure (maybe)
% +++++ colors in Raster ++++
for i=1:length(labels_frames)
    AC=Experiment(signif_frames(i),indexes); % Active Cells
    AC=find(AC==1);
    for j=1:length(AC)
        % k=indexes(j);
        subplot(3,1,[1,2]); hold on;  % ***************Raster
        plot((1/(fs))*signif_frames(i)/60,AC(j),'.',...
                'Color',ColorState(labels_frames(i),:),...
                'MarkerSize',10); 
%             'MarkerEdgeColor',ColorState(labels_frames(i),:),...
%             'MarkerFaceColor',ColorState(labels_frames(i),:),...
    end
    subplot(3,1,3); hold on;
    Cframe=signif_frames(i);
    if i+1>length(labels_frames) % for the last frame
        Nframe=Cframe;
    else
        Nframe=signif_frames(i+1);
    end
    if Nframe-Cframe>1
        plot((1/(fs))*(Cframe-1)/60,CoAc(Cframe),...
            '.','Color',ColorState(labels_frames(i),:),...
            'MarkerSize',8)
    elseif Nframe-Cframe==1
        plot((1/(fs))*[Cframe-1,Nframe-1]/60,[CoAc(Cframe),CoAc(Nframe)],...
            '.','Color',ColorState(labels_frames(i),:),...
            'LineWidth',5)
    end

end
% Coactivity *******************************************************

subplot(3,1,3); hold on; % ******************** Coactivity Threshold
plot((1/fs)*[0,length(CoAc)]/60,[THR,THR],'--','Color','k')
ax=gca;
MaxYY=max(sum(Experiment,2));
if MaxYY<=THR
    MaxYY=THR+1;
end
ax.YTick=[THR,MaxYY];
axis([0,(1/fs)*(length(CoAc))/60,0,max(sum(Experiment,2))+1])