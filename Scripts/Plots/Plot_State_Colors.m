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
        rectangle('Position',[xposition(1),ypositon(1),...
                xposition(2),ypositon(2)],'Curvature',0.2,...
                'EdgeColor',ColorState(labels_frames(i),:),...
                'FaceColor',ColorState(labels_frames(i),:));
        fprintf('>>coloring cell %i ',AC(j))
    end
%     drawnow; (slows stuff)
    fprintf('\n')
end
%% CAG  *******************************************************
disp('>>Coloring Ensembles at CAG')
subplot(3,1,3); hold on; % ******************** Coactivity Threshold
f=1;
while f<=numel(signif_frames)
    % Get Segment of the Same Color
    if f<numel(signif_frames)
        LengthColor=f+1;
    else
        LengthColor=f;
    end
    frameA=f;
    frameB=frameA;
    if or(LengthColor>numel(signif_frames),f>numel(signif_frames))
        frameB=0;
    else
        aux1=0;
        while labels_frames(f)==labels_frames(LengthColor) && ...
                signif_frames(f)==signif_frames(LengthColor)-1 &&...
                f<=numel(labels_frames) && LengthColor<=numel(labels_frames)
            frameB=LengthColor-frameA;
            f=f+1;
            aux1=aux1+1;
            fprintf('*')
        end
        if aux1==0
            frameB=0;
        end
    end
    
    % Plot Segment
    %Segment=[frameA,frameA+frameB];
    if frameB>0
        % LINE
        plot(signif_frames(frameA:frameA+frameB)/fs/60,CoAc(signif_frames(frameA:frameA+frameB)),...
        'Color',ColorState(labels_frames(frameA),:),...
        'LineWidth',2);
    else
        % DOT
        plot(signif_frames(frameA)/fs/60,CoAc(signif_frames(frameA)),...
            'Marker','.','MarkerSize',3.5,...
            'MarkerEdgeColor',ColorState(labels_frames(frameA),:),...   
            'MarkerFaceColor',ColorState(labels_frames(frameA),:));
    end
    
    f=f+1;
    fprintf('\n')
end
disp('>>Colors at CAG.')