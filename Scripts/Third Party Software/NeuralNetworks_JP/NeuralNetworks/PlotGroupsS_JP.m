% Plot groups spatially
% Plot the different groups spatially.
%
% PlotGroupsS(XY, Active, CellsGroups, numFig)
%
% Inputs
%
%
% Output
% 
%
% ..:: by Jesús E. Pérez-Ortega ::.. Jun-2012
%
% V 2.0 Replacement "g=length(g2plot)" by this "g=max(g2plot)" jun-2012

function PlotGroupsS_JP(XY, Active, CellsGroups, g2plot, numFig)

g=max(g2plot);
ptsymb = {'.r','.g','.b','.c','.m','.y','.k'};

figure(numFig)

for i=1:g
    if find(g2plot==i)
        subplot(1,g,i)
        plot(XY(:,1),XY(:,2),'O','color',[.5 .5 .5])
        view([0 -90])
        hold on
        
        cells=find(CellsGroups(i,:));
        ActiveCells=Active(cells);
        plot(XY(ActiveCells,1),XY(ActiveCells,2),ptsymb{i},'MarkerSize',18)
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        hold off
    end
end


