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
function Plot_State_Colors(labels_frames,signif_frames,ColorState,Experiment,fs,CoAc,indexes)
% Search of figure (maybe)
TransientHight=0.75; % Same as Plot_Raster_Ensembles.m
% +++++ colors in Raster +++++
for i=1:length(labels_frames)
    AC=Experiment(signif_frames(i),indexes); % Active Cells
    AC=find(AC==1);
    for j=1:numel(AC)
        % k=indexes(j);
        subplot(3,1,[1,2]); hold on;  % ***************Raster
%         plot((1/(fs))*signif_frames(i)/60,AC(j),'Marker','square',...
%                 'MarkerEdgeColor',ColorState(labels_frames(i),:),...
%                 'MarkerFaceColor',ColorState(labels_frames(i),:),...
%                 'MarkerSize',3.4);
        xposition=(1/(fs))*[signif_frames(i)-0.5,1]/60;
        ypositon=[AC(j)-TransientHight/2,TransientHight];
        rectangle('Position',[xposition(1),ypositon(1),...
                xposition(2),ypositon(2)],'Curvature',0.2,...
                'EdgeColor',ColorState(labels_frames(i),:),...
                'FaceColor',ColorState(labels_frames(i),:));
        fprintf('°')
    end
    fprintf('\n')
%     % CAG COLORS **********************
%     subplot(3,1,3); hold on;
%     Cframe=signif_frames(i);
%     if i+1>length(labels_frames) % for the last frame
%         Nframe=Cframe;
%     else
%         Nframe=signif_frames(i+1);
%     end
%     if Nframe-Cframe>1
%         plot((1/(fs))*(Cframe-1)/60,CoAc(Cframe),...
%             'Color',ColorState(labels_frames(i),:),...
%             'LineWidth',2)
%     elseif Nframe-Cframe==1
%         plot((1/(fs))*[Cframe-1,Nframe-1]/60,[CoAc(Cframe),CoAc(Nframe)],...
%             'Marker','.',...
%             'MarkerEdgeColor',ColorState(labels_frames(i),:),...
%             'MarkerFaceColor',ColorState(labels_frames(i),:),...
%             'MarkerSize',3);
%     end

end


% CAG  *******************************************************
disp('>>Coloring Ensmebles at CAG')
% subplot(3,1,3); hold on; % ******************** Coactivity Threshold
% plot((1/fs)*[0,length(CoAc)]/60,[THR,THR],'--','Color','k')
% ax=gca;
% MaxYY=max(sum(Experiment,2));
% if MaxYY<=THR
%     MaxYY=THR+1;
% end
% ax.YTick=[THR,MaxYY];
% axis([0,(1/fs)*(length(CoAc))/60,0,max(sum(Experiment,2))+1])
disp('>>Colors at CAG.')