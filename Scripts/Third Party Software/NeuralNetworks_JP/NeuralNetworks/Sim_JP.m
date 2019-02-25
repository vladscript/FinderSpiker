% Similarity of peaks
%
% Get the similarity between vectors of peaks.
%
% [Sim] = Sim_JP(X, method, numFig)
%
% Inputs
% X = binary data as matrix PxC (P = #peaks, C = #cells)
% method = string of the method to use ('cosine, 'jaccard', 'euclidean') 
% numFig = number of the figure to plot
% 
% Outputs
% Sim = similarity as matrix PxP
%
% ..:: by Jesús E. Pérez-Ortega ::.. Jun-2012

function [Sim] = Sim_JP(X, method, numFig)

dist=pdist(X,method);
Sim=1-squareform(dist);

figure(numFig)
pcolor(Sim)
title('Similarity')
end
