% Hierarchical Tree
% Get hirarchical tree from de similarity index of vectors.
%
% [Tree] = HTree_JP(Sim, numFig)
%
% Inputs
% Sim = similarity as matrix PxP (P=#peaks)
% numFig = number of the figure to plot
% method = method to create a hierarchical cluster ('average','centroid',
%           'complete','median','single','ward','weighted')
%
% Output
% Tree = hirarchical tree.
%
% ..:: by Jesús E. Pérez-Ortega ::.. Jun-2012
%
% V 2.0 3rd input added: 'method' Aug-2012

function [Tree] = HTree_JP(Sim, numFig, method)

% Get Dissimilarity (or distance)
dis=1-Sim;
peaks=size(Sim,1);

% Hierarchical Cluster Tree
Tree = linkage(squareform(dis,'tovector'),method);

% Plot Similarity and Hierarchical Tree
figure(numFig)

subplot(1,4,4)
dendrogram(Tree,0,'orientation','rigth');
view([0 90])
Tid=str2num(get(gca,'YTicklabel'));
subplot(1,4,1:3)
pcolor(Sim(Tid,Tid))
set(gca,'YTicklabel',Tid)
set(gca,'XTicklabel',Tid)
set(gca,'YTick',[1:peaks])
set(gca,'XTick',[1:peaks])
title(method)