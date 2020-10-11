% Plot Raster Sorted by Activity Rate:
%  Retieve Variable: >>IndexesbySubgroups, sorted by:
% Inhibited
% Depressed
% Unchanged
% Facilitated
% Activated
% Input
%   R_Condition:     Cell of ACtivity Matrices
%   Condition_Names: Cell of Condition Strings
%   fs:              Sampling Frequency
%   Condition to Compare: selected within
% Output
%   NewIndex:        Sorted by Activity Diffrence
function NewIndex=plot_activitychange(R_Condition,Condition_Names,fs)

N=numel(R_Condition);
R_Total=[];
for n=1:N
    R=R_Condition{n};
    R_Total=[R_Total,R];
    AR(:,n)=sum(R,2)/size(R,2);
end
Ncells=size(R,1);

[ActivityGroups, ~, ~, RoAG]= vennRASTER(R_Condition,Condition_Names);

%% RE-SORT
% Depressed
[~,IndxDep]=sort(RoAG.dep(:,2));
% Unchanged (a 10 %)
[~,IndxUn]=sort(RoAG.unc(:,2)-RoAG.unc(:,1));
% Facilitated
[~,IndxFac]=sort(RoAG.fac(:,1));

NewIndex=[ActivityGroups.Ndep(IndxDep);ActivityGroups.Nunc(IndxUn); ...
    ActivityGroups.Nfac(IndxFac(end:-1:1))];

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

