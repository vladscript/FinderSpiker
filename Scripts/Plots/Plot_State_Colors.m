%% Plotting States Colors  *********************************************
% Requeriment:
% Plot Raster
% Input
% 
%   labels_frames: Ensemble Labels in each frame
%   signif_frames: Significative Coactivity Frames
%   ColorState:    Color of The Ensembles
%   Experiment:    Raster Concatenated and Ensemble Sorted
%   fs:            Sampling Frequency
%   CoAc:          Coactivity signal
%   indexes:       Sorting
% 
% Ouput
% 
%   Plot on the Raster plot
function Plot_State_Colors(labels_frames,signif_frames,ColorState,Experiment,fs,CoAc,indexes)
%% Search of figure (maybe)
TransientHight=0.75; % Same as Plot_Raster_Ensembles.m
%% +++++ colors in Raster +++++
subplot(3,1,[1,2]); hold on;  % ***************@ Actual Raster
for i=1:length(labels_frames)
    AC=Experiment(signif_frames(i),indexes); % Active Cells
    AC=find(AC==1);
    for j=1:numel(AC)
        % k=indexes(j);
        xposition=(1/(fs))*[signif_frames(i)-0.5,1]/60;
        ypositon=[AC(j)-TransientHight/2,TransientHight];
        if labels_frames(i)>0
            rectangle('Position',[xposition(1),ypositon(1),...
                    xposition(2),ypositon(2)],'Curvature',[0,0],...
                    'EdgeColor',ColorState(labels_frames(i),:),...
                    'FaceColor',ColorState(labels_frames(i),:));
        end
        fprintf('>>coloring cell %i ',AC(j))
    end
%     drawnow; (slows stuff)
    fprintf('\n')
end

%% CAG  *******************************************************
% 
% USELESS CODE: color intervals are plotted according to the smoothed
% version of the CAG signal when Ensemble intervals are computed.
% 
% disp('>>Coloring Ensembles at CAG')
% subplot(3,1,3); hold on; % ******************** Coactivity Plot
% f=1; % additional consecutive label frames
% while f<=numel(signif_frames)
%     % Get Segment of the Same Color
%     if f<numel(signif_frames)
%         frameB=f+1;
%     else
%         frameB=f; % Limit
%     end
%     frameA=f;
%     % Get consecutive frames of the same ensemble
%     reachlimint=true;
%     pref=f;
%     auxinloop=1;
%     while reachlimint;
%         if labels_frames(frameA)==labels_frames(frameB) && ... % ensemble
%             signif_frames(frameA)==signif_frames(frameB)-auxinloop     % consecutive frame
%             f=f+1;
%             auxinloop=auxinloop+1;
%             fprintf('~')
%         else
%             reachlimint=false;
%             if auxinloop>1
%                 f=f-1;
%             end
%         end
%         if frameB>numel(signif_frames)
%             reachlimint=false;
%         end
%         frameB=f;
%     end
%     
%     if pref==f
%         frameB=frameA;
%     end
%     % Plot Segment
%     if labels_frames(frameA)>0
%         if frameB>frameA
%             plot(signif_frames(frameA:frameB)/fs/60,CoAc(signif_frames(frameA:frameB)),...
%             'Color',ColorState(labels_frames(frameA),:),...
%             'LineWidth',1);
%         elseif frameB==frameA
%             DotA=signif_frames(frameA)-0.5;
%             DotB=signif_frames(frameA)+0.5;
%             if DotB>numel(signif_frames)
%                 DotB=numel(signif_frames);
%             end
%             CAG_A=signif_frames(frameA)-1;
%             CAG_B=signif_frames(frameA)+1;
%             if CAG_A<1
%                 CAG_A=1;
%             end
%             if CAG_B>numel(signif_frames)
%                 CAG_B=numel(signif_frames);
%             end
%             plot([DotA/fs/60,DotB/fs/60],...
%                 [mean([CoAc(CAG_A),CoAc(signif_frames(frameA))]),...
%                 mean([CoAc(signif_frames(frameA)),CoAc(CAG_B)])],...
%             'Color',ColorState(labels_frames(frameA),:),...
%             'LineWidth',1);
%         end
%         % drawnow;
%     end
%     
%     f=f+1;
%     fprintf('\n')
% end
% disp('>>Colors at CAG.')