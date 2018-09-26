% Higlight Neuron in Raster Plot in Red
% Input
%   Plot of Raster:     Current Figure (from Raster_Plot_V)
%   Neuron ID           Neuron ID
% Ouput
% @ raster plot:
%       Rectangle
function HighlightNeuronatRaster(NeuronID)
% Example:
axisfig=gcf;
NeuronIDstr=num2str(NeuronID);
LengthStr=length(NeuronIDstr);
% Do:
TicksY=axisfig.Children(2).YTick;           % LOCATIONS
LabelsY=axisfig.Children(2).YTickLabel;     % cell IF ACTUAL SORTING
NewTick=NeuronID; % Tick-index of Neuron
LabelsNum=str2num(LabelsY);

if ~ismember(NeuronID,LabelsNum) % check it's NOT already in Yticks        
    TicksY=sort([TicksY,NewTick]);
    index4label=find(TicksY==NewTick);
    if index4label<length(TicksY) % if it's not the last
        movingLabels=LabelsNum(index4label:end)';
        %         LabelsY(index4label)=AddIndx;
        LabelsNum(index4label)=NeuronID;
        LabelsNum=[LabelsNum(1:index4label);movingLabels'];
    else
        LabelsNum=[LabelsNum';NeuronID];
    end
else
    %LabelsNum=LabelsNum';
    disp('Already in PLOT')
end
TickLocation=find(LabelsNum==NewTick);
axisfig.Children(2).YTickLabel=num2str( LabelsNum );
% Update Ticks and Label Ticks
axisfig.Children(2).YTick       =   TicksY;
% axisfig.Children(2).YTickLabel  =   LabelsY;
% Xlims of subplot
Xlim=axisfig.Children(2).XLim;
% Highlight Rectangle
subplot(3,1,[1,2]); hold on; 
rectangle('Position',[Xlim(1),LabelsNum(TickLocation)-0.25,Xlim(2),0.5],...
    'EdgeColor',[0.5,0.5,0.5],...
    'LineStyle','--',...
    'LineWidth',1.15);
hold off;
