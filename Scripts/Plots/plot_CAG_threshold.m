% Use Only if there is a Raster Plotted Already
% Input:
%   THR: Cell of Threshold(s) by Condition
%   Cell of Raster Lengths (frames) by Condition
%   Sampling Frequency (Hz)
% Outpu:
% Segmented line of THreshold in subplot(1,3,3)
function plot_CAG_threshold(THR,LENGHTRASTER,fs)
[~,NC]=size(LENGHTRASTER);
ax=subplot(3,1,3); 
hold(ax,'on'); 
ax.YTick=[];
% Read & Get Proportion of Raster Lengths
FramesTotal=0;
FramesConditionCumm=0;
% Experiment=[];
for i=1:NC
    FramesCondition=LENGHTRASTER{i};
    % Experiment=[Experiment,LENGHTRASTER{i}];
    FramesConditionCumm=FramesConditionCumm+FramesCondition;
    plot((1/fs)*[FramesTotal,FramesConditionCumm]/60,[THR{i},THR{i}],'--','Color','k')
    % ax=gca;
    MaxYY=max(ax.Children(end).YData);
    if MaxYY<=THR{i}
        MaxYY=THR{i}+1;
    end
    ax.YTick=sort(unique([ax.YTick,THR{i},MaxYY]));
    FramesTotal=FramesCondition;
end
% axis([0,(1/fs)*(FramesConditionCumm)/60,0,max(sum(Experiment))+1])