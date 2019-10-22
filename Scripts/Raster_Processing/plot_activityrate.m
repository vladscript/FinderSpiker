% Plot Raster Sorted by Activity Rate:
% AR=Active Frames/Total Frames by Condition
% Run after Selected Conditions
% Input
%   R_Condition:     Cell of ACtivity Matrices
%   Condition_Names: Cell of Condition Strings
%   fs:              Sampling Frequency
%   Condition to Compare: selected within
% Output
%   NewIndex:        Sorted by Activity Diffrence
function NewIndex=plot_activityrate(R_Condition,Condition_Names,fs)
N=numel(R_Condition);
R_Total=[];
for n=1:N
    R=R_Condition{n};
    R_Total=[R_Total,R];
    AR(:,n)=sum(R,2)/size(R,2);
end
Ncells=size(R,1);
% Create Figure For the Experiment
Plot_Raster_Ensembles(R_Total,fs);
Label_Condition_Raster(Condition_Names,R_Condition,fs);   % Labels
RasterFig=gcf;
RasterFig.Children(end)
ARfig=figure;
ARfig.Position=[403 246 200 420];
ax1=subplot(3,1,[1,2]);
plot(ax1,AR,repmat([1:Ncells]',1,N),'LineWidth',2)
ax1.YLim=RasterFig.Children(end).YLim;
ax1.YTick=RasterFig.Children(end).YTick;
legend(ax1,Condition_Names,'Location','southoutside');
linkaxes([ax1,RasterFig.Children(end)],'y')
fprintf('>>\n\n You can move the legend box using the cursor\n\n')

% Sorting
% Choose Experimental Conditon: 2,3,...k,
OrdinalString={'1st','2nd'};
if N>1
    for c=1:2
        [index_var{c},index_CHECK{c}] = listdlg('PromptString',['Select ',OrdinalString{c},' Condition'],...
                        'SelectionMode','single',...
                        'Name',['>> Condition: [',num2str(c),'/2] '],...
                        'ListSize',[300 200],...
                        'ListString',Condition_Names);            
    end
    fprintf('>>\nComparing: %s  vs %s\n',Condition_Names{index_var{1}},Condition_Names{index_var{2}})
    % [ Condition 2- Condition 1 ] (Activity Rate)
    ARdiff=diff(AR(:,[index_var{1},index_var{2}]),1,2);
    % Define Two Groups: Augmented-Diminished
    AugmentedSet=find(ARdiff>=0);
    DiminishedSet=find(ARdiff<0);
    % Define Two SubGroups: 
    %     Augmented:{Less More} 
    % [~,AugSortedIndx]=sort(ARdiff(AugmentedSet));
    [~,AugSortedIndx]=sort(AR(AugmentedSet,index_var{2}));
    
    %     Diminished {Less to More]
    % [~,DimSortedIndx]=sort(ARdiff(DiminishedSet));
    [~,DimSortedIndx]=sort(AR(DiminishedSet,index_var{1}),'descend');
    
    NewIndex=[AugmentedSet(AugSortedIndx);DiminishedSet(DimSortedIndx)];
else % Only Sort by AR
    [~,NewIndex]=sort(AR);
end
% PLot Sorted Raster
Plot_Raster_Ensembles(R_Total,fs,1,NewIndex);
Label_Condition_Raster(Condition_Names,R_Condition,fs);   % Labels
%  Replot AR
RasterFig2=gcf;
RasterFig2.Children(end);
ARfig=figure;
ARfig.Position=[503 246 200 420];
ax2=subplot(3,1,[1,2]);
plot(ax2,AR(NewIndex,:),repmat([1:Ncells]',1,N),'LineWidth',2)
ax2.YLim=RasterFig2.Children(end).YLim;
ax2.YTick=RasterFig2.Children(end).YTick;
ax2.YTickLabel=RasterFig2.Children(end).YTickLabel;
ax2.XLabel.String='ActiveFrames/TotalFrames';
ax2.XLabel.FontSize=9;
legend(ax2,Condition_Names,'Location','southoutside');
linkaxes([ax2,RasterFig2.Children(end)],'y');
fprintf('>>\n\n You can move the legend box using the cursor\n\n')

