%% Set Labels to Whole Raster: 
% Run Only After Plot_Raster_V(RASTER_WHOLE_Clean,fs);
% Input
%   Cell of Names of Conditons
%   Raster Lengths
%   Sampling Frequency
% Output
%    Annotation to Current Figure in Minutes
function Label_Condition_Raster(Names_Conditions,Raster_Condition,fs)
[~,NC]=size(Raster_Condition);

% Read & Get Proportion of Raster Lengths
FramesTotal=0;
for i=1:NC
    [~,FramesCondition(i)]=size(Raster_Condition{i});
    FramesTotal=FramesTotal+FramesCondition(i);
end
FractionFrame=FramesCondition/FramesTotal;
% Position Text Limits in X-axis
InitialX_Text=0.13;
FinalX_Text=0.93;
DeltaX=FinalX_Text-InitialX_Text;
% Widths of Box Texts by a 3-rule: DeltaX->FramesTotal
Widthboxtex=FractionFrame*DeltaX;
% Plotting Stuff
SumFrames=1;
figure1=gcf;
Xboxtex=InitialX_Text;
Yboxtex=0.93;
Heightboxtex=0.04;
for i=1:NC
    [Cells,~]=size(Raster_Condition{i});    
    Frames=FramesCondition(i);
    SumFrames=SumFrames+Frames;
    Position_Minutes=SumFrames/60/fs;
    % Create textbox
    % Position x,y,w,h
    annotation(figure1,'textbox',...
        [Xboxtex Yboxtex Widthboxtex(i) Heightboxtex],...
        'String',Names_Conditions{i},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',12,...
        'FontName','Arial',...
        'FitBoxToText','off');
    Xboxtex=Xboxtex+Widthboxtex(i);
    % Create line
    subplot(3,1,[1,2]); hold on;
    line([Position_Minutes,Position_Minutes],[-Cells,Cells],...
        'Color',[0.8 0.8 0.8],...
        'LineWidth',2,...
        'LineStyle',':');
    subplot(3,1,3); hold on;
    line([Position_Minutes,Position_Minutes],[-Cells,Cells],...
        'Color',[0.8 0.8 0.8],...
        'LineWidth',2,...
        'LineStyle',':');
end