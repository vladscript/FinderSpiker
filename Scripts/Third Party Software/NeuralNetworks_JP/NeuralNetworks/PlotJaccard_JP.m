% Jaccard correlation
% Compute and plot the jaccard correlation between cells.
%
% [Jaccard AdjacencyMatrix] = PlotJaccard_JP(X, XY, th, sameWidth, colorCells, numFig)
%
% Inputs
% X = data as F x C matrix (F = #frames, C = #cells)
% XY = data as C x 2 matrix (1st & 2nd columns are X&Y coordinates
% th = correlation index threshold for drawing significant edges
% sameWidth = indicate if you want to draw all edges with the same width or
%             with different width depending of correlation index (1 =
%             same; 0 = different)
% colorCells = color of cells as vector 1x3 [R G B]
% numFig = number of the figure to plot
%
% Output
%
% Jaccard = Jaccard correlation index between all cells as C x C matrix
% AdjacencyMatrix = Adjacency matrix of significant correlations
%
% ..:: by Jesús E. Pérez-Ortega ::.. Jun-2012
%
% V 2.0 5th input added: 'sameWidth'
%       6th input added: 'colorCells' Jun-2012

function [Jaccard AdjacencyMatrix]=PlotJaccard_JP(X, XY, th, sameWidth, colorCells, numFig)

Jaccard=1-squareform(pdist(X','jaccard'));
AdjacencyMatrix=Jaccard>th;

C=size(X,2);
figure(numFig)
hold on
for a=1:C-1
    for b=a+1:C
        if Jaccard(a,b)>th
            if sameWidth
                plot(XY([a b],1),XY([a b],2),'color',[.5 .5 .5],'LineWidth',1);
            else
                plot(XY([a b],1),XY([a b],2),'color',[.5 .5 .5],'LineWidth',Jaccard(a,b)*20);
            end
        end
    end
end
plot(XY(:,1),XY(:,2),'.','color',colorCells,'MarkerSize',23)
title('Significant correlation')
view([0 -90])
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])