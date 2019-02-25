% Jaccard Montecarlo
% Get the probability distribution of cells correlated from raster
% by Montecarlo method.
%
% [th P bins] = MCJaccard_JP(X, i, numFig)
%
% Inputs
% X = data as F x C matrix (F = #frames, C = #cells).
% i = number of simulations.
% alpha = probability of error.
% numFig = number of the figure to plot.
%
% Output
% th = threshold from probability distribution < alpha
% P = probability distribution of correlation index
% bins = probability distribution intervals
%
% ..:: by Jesús E. Pérez-Ortega ::.. jun-2012

function [th P bins] = MCJaccard_JP(X,i,alpha,numFig)

% Get size of data
[F C] = size(X);

% Get jaccard correlation index
Jaccard=1-pdist(X','jaccard');
n=length(Jaccard);
bin=100;
bins=0:1/bin:1;

% Observed data
HO = histc(Jaccard,bins);   % Histogram of correlation index
PO = HO/n;                      % Correlation probability (OBSERVED)

% Simulated data
PS=zeros(i,bin+1);
for j = 1:i                             % Make i simulations
    for k = 1:C                           
        XS(:,k) = X(randperm(F),k);     % Shuffle spikes for every cell
    end
  
    JaccardSim=1-pdist(XS','jaccard');
    HS = histc(JaccardSim,bins);        % Histogram of correlation index (simulated)
    PS(j,:) = HS/n;                     % Correlation probability (simulated)
end
PE=mean(PS);                            % Mean of Probability (EXPECTED)

% Get propability
P=zeros(1,bin+1);
for j = 1:i
    P = P+double(PO>=PS(j,:));
end
P = 1-P/i;

% Set threshold from Probability distribution < alpha
th=find(P(4:101)<alpha,1,'first')+3;
interval=1:th+5;
th=bins(th);
% Probability of get Observed data
figure(numFig)
plot(bins(interval),PO(interval),'o-r')
title('Probability of obtaining observed data')
hold on
plot(bins(interval),PE(interval),'o-b')
plot(bins(interval),PO(interval)-PE(interval),'o-g')
plot(bins(interval),P(interval)','o-M')
hold off
legend('Observed', strcat('Expected (',num2str(i),' simulations)'),...
    'Substract ( Observed - Expected )','Probability')