% Get Clusters
% Get clusters from the hierarchical cluster tree
%
% [PeaksGroups CellsGroups] = Cluster_JP(X, Tree, Groups)
%
% Inputs
% X = binary data as matrix PxC (P = #peaks, C = #cells)
% SimOrXpeaks = Sim, similarity as matrix PxP (P=#peaks) for metrics Dunn &
%               Contrast; Xpeaks, peaks vectors as matrix PxC for metrics
%               Conectivity & Davies
%               (P = #peaks; C=#cells)
% Metric = index to compute ('Dunn','Conectivity','Davies','Contrast')
% numFig = number of the figure to plot
%
% Output
% ClustIdx = clustering indexes of 'Metric' from 2 to 10 groups
% Groups = best number of groups according to the index selected
%
% ..:: by Jesús E. Pérez-Ortega ::.. Jun-2012

function [PeaksGroups CellsGroups] = Cluster_JP(X, Tree, Groups)

% Plot Raster with only peaks per group sorted
[Fp Cp] = size(X);
PeaksGroups = cluster(Tree,'maxclust',Groups);

% Description of groups
CellsGroups=zeros(Groups, Cp);
for i=1:Groups
    PeaksGroupIdx=find(PeaksGroups==i);
    PeaksGroup=numel(PeaksGroupIdx);
    CellsGroup=find(sum(X(PeaksGroupIdx,:),1));
    CellsGroups(i,CellsGroup)=1;
    disp(['     Group ' num2str(i) ': ' num2str(numel(CellsGroup)) ' cells'...
        ' (' num2str(PeaksGroup) ' peaks) (cells: ' num2str(CellsGroup)])
end
disp(['     Total cells: ' num2str(Cp)])