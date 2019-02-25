% Clustering indexes
% Get indexes for evaluating clustering from hierarchical cluster tree
%
% [ClustIdx g] = ClustIdx_JP(Tree, SimOrXpeaks, Metric, numFig)
%
% Inputs
% Tree = hierarchical cluster tree
% SimOrXpeaks = Sim, similarity as matrix PxP (P=#peaks) for metrics Dunn &
%               Contrast; Xpeaks, peaks vectors as matrix PxC for metrics
%               Connectivity & Davies
%               (P = #peaks; C=#cells)
% Metric = index to compute ('Dunn','Connectivity','Davies','Contrast')
% numFig = number of the figure to plot
%
% Output
% ClustIdx = clustering indexes of 'Metric' from 2 to 10 groups
% g = best number of groups according to the index selected
%
% ..:: by Jesús E. Pérez-Ortega ::.. Jun-2012

function [ClustIdx g] = ClustIdx_JP(Tree, SimOrXpeaks, Metric, numFig)

Dist=1-SimOrXpeaks;

for i=2:10
    T = cluster(Tree,'maxclust',i);
    idxSort=[];
    g=max(T);
    for i=1:g
        Tidx2=find(T==i);
        idxSort=[Tidx2; idxSort];
    end
    
    switch(Metric)
        case 'Dunn'
            ClustIdx(g)=DunnIdx_JP(g,Dist,T);
        case 'Connectivity'
            ClustIdx(g)=CCidx2_JP(g,SimOrXpeaks,T);
        case 'Davies'
            ClustIdx(g)=DaviesIdx_JP(g,SimOrXpeaks,T);
        case 'Contrast'
            ClustIdx(g)=ContrastIdx_JP(g,SimOrXpeaks,T);
    end
end
figure(numFig)
plot(ClustIdx)
hold on

if strcmp(Metric,'Contrast')
   [gIdx g]=max(ClustIdx);
   gIdx=gIdx*.9;
   g=find(ClustIdx>gIdx,1,'first');
   gIdx=ClustIdx(g);
else
   [gIdx g]=max(ClustIdx);
end

plot(g,gIdx,'*r')
hold off
title([Metric '''s index (' num2str(g) ' groups recommended)'])
