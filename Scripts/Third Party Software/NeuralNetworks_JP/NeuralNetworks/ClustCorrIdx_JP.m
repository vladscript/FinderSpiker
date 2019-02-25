% Clustering Correlation Indexes
% Get teh correlation indexes for evaluating clustering
% from hierarchical cluster tree.
%
% [ClustIdx g] = ClustCorrIdx_JP(Tree, CorrelationMatrix, numFig)
%
% Inputs
% Tree = hierarchical cluster tree
% CorrelationMatrix = correlation as a square matrix
% numFig = number of the figure to plot
%
% Output
% ClustIdx = clustering indexes from 2 to 10 groups
% g = best number of groups according to the index selected
%
% ..:: by Jesús E. Pérez-Ortega ::.. Aug-2012
% V2.0 Replacement idx_group "(...*n-1)/(n-1)" by "()/n" modified Nov-2012

function [ClustIdx g] = ClustCorrIdx_JP(Tree, Data, CorrelationMatrix, numFig)
CorrelationMatrix=CorrelationMatrix+eye(size(CorrelationMatrix));
ClustIdx=zeros(10,1);
for i=2:10
    T = cluster(Tree,'maxclust',i);
    g=max(T);
    
    idx_group=zeros(1,g);
    for j=1:g
        group=T==j;
        idx=find(sum(Data(group,:),1));
        n=length(idx);
        idx_group(j)=(mean(mean(CorrelationMatrix(idx,idx))))/n;
    end
    ClustIdx(i)=mean(idx_group)/std(idx_group);
    g=max(ClustIdx);
end

figure(numFig)
plot(ClustIdx,'-b');
hold on
[gIdx g]=max(ClustIdx);
plot(g,gIdx,'*r')
title(['Correlation''s index (' num2str(g) ' groups recommended)'])