%% Script to fit linear accumulated growing rate of activity for each raster
% Input
%   Cell of Raster for each Condition
% Output
%   M: vector of slope for each condition
%   B: vector of y-intercept for each condition
% display summarize table
% example:
% >> [M,B]=slope_that_raster(R_Condition,varargin)
function [M,B]=slopes_from_raster(R_Condition,fs,Names_Conditions,varargin)
%% setup
NC=numel(R_Condition); % Number of Conditons
% Color Selector
[CM,ColorIndx]=Color_Selector(Names_Conditions);
if numel(varargin)>0
    MetaDataColocaliation=varargin{1};
    fprintf('\n>>Analyze for %s + cells detected with %s\n',...
        MetaDataColocaliation.Cells{1},MetaDataColocaliation.Dye{1});
    CellKind=MetaDataColocaliation.Cells{1}(isstrprop(MetaDataColocaliation.Cells{1}, 'alphanum'));
    PositiveCells=MetaDataColocaliation.PositiveCells;
    NegativeCells=MetaDataColocaliation.NegativeCells;
    AreMerge=true;
    LineObjectsPos={};
    LineObjectsNeg={};
else
    fprintf('\n>>Without specific cell type\n');
    AreMerge=false;
    CellKind='All';
end
%% Create figure objects ***************************************
figAll=figure('NumberTitle','off','Name','All Cells');
if AreMerge
    figPos=figure('NumberTitle','off','Name',['+ ',CellKind,' Cells']);
    figNeg=figure('NumberTitle','off','Name',['- ',CellKind,' Cells']);
end
%% MAIN LOOP
SLOPE=zeros(NC,3); % size: Nconditions vs [ALL/+/-]
Yinter=zeros(NC,3);
LineObjects={};
for n=1:NC
    % ALL ***************************************************
    R=R_Condition{n};
    p=slopethatraster(R,fs);
    % PLOT *************************************************    
    figure(figAll)
    h=plot(1:size(R,2),cumsum(sum(R)),'LineWidth',3,...
    'color',CM(ColorIndx(n),:)); hold on
    plot(1:size(R,2),p(1)*[0:size(R,2)-1]/fs+p(2),'LineWidth',1,...
    'color','k'); hold on
    axis tight; grid on;
    LineObjects=[LineObjects;h];
    if AreMerge
        % POSITIVE ******************************************
        Rpos=R_Condition{n}(PositiveCells,:);
        pp=slopethatraster(Rpos,fs);
        % NEGATIVE ******************************************
        Rneg=R_Condition{n}(NegativeCells,:);
        pn=slopethatraster(Rneg,fs);
        % PLOTS *************************************************    
        figure(figPos)
        hpos=plot(1:size(Rpos,2),cumsum(sum(Rpos)),'LineWidth',3,...
        'color',CM(ColorIndx(n),:)); hold on
        plot(1:size(Rpos,2),pp(1)*[0:size(Rpos,2)-1]/fs+pp(2),'LineWidth',1,...
        'color','k'); hold on
        axis tight; grid on;
        LineObjectsPos=[LineObjectsPos;hpos];
        figure(figNeg)
        hneg=plot(1:size(Rneg,2),cumsum(sum(Rneg)),'LineWidth',3,...
        'color',CM(ColorIndx(n),:)); hold on
        plot(1:size(Rneg,2),pn(1)*[0:size(Rneg,2)-1]/fs+pn(2),'LineWidth',1,...
        'color','k'); hold on
        axis tight; grid on;
        LineObjectsNeg=[LineObjectsNeg;hneg];
    else
        pp=0; pn=0;
    end
    % MAKE TABLE ##################################################
   SLOPE(:,n)=[p(1);pp(1);pn(1)];
   Yinter(:,n)=[p(end);pp(end);pn(end)];
   NamesOK{n}=Names_Conditions{n}(isstrprop(Names_Conditions{n}, 'alphanum'));
end
%% Fix axis
figure(figAll)
AxisCumSum=gca;
fixaxis(AxisCumSum,fs)
legend(LineObjects,NamesOK)
if AreMerge
    figure(figPos)
    AxisCumSum=gca;
    fixaxis(AxisCumSum,fs)
    legend(LineObjectsPos,NamesOK)
    figure(figNeg)
    AxisCumSum=gca;
    fixaxis(AxisCumSum,fs)
    legend(LineObjectsNeg,NamesOK)
end
% Display OUTPUT:
M=mat2dataset(SLOPE,...
    'ObsNames',{'All',['Positive_',CellKind],['Negative_',CellKind]},...
    'VarNames',NamesOK)
B=mat2dataset(Yinter,...
    'ObsNames',{'All',['Positive_',CellKind],['Negative_',CellKind]},...
    'VarNames',NamesOK)
% NESTED functions
function fixaxis(AxisCumSum,fs)
    AxisCumSum.YLabel.String='Cummulative CAG sum [a.u.]';
    AxisCumSum.XLabel.String='Minutes';
    Xlimit=AxisCumSum.XLim(end);
    NewXlTicks=[0:60*fs:Xlimit];
    AxisCumSum.XTick=NewXlTicks;
    for naxis=1:numel(NewXlTicks)
        AxisCumSum.XTickLabel{naxis}=num2str(NewXlTicks(naxis)/60/fs);
    end
end
end