% Raster Montecarlo 3 (cummulative)
% Get the probability distribution of cells sincronized from raster
% by Montecarlo method.
%
% [th P] = RMC(X, i, alpha, numFig)
%
% Inputs
% X = data as F x C matrix (F = #frames, C = #cells).
% i = number of simulations.
% alpha = probability of error.
% numFig = number of the figure to plot.
%
% Output
% P = probability distribution along C as row vector
% th = threshold from probability distribution < alpha
%
% ..:: by Jesús E. Pérez-Ortega ::.. Oct-2011
%
% V 2.0 Replacement "interval=[0:th]" by this "interval=[0:th+5]" Jun-2012
% V 3.0 Cummulative probability, asterisks in significative points


function [th P] = MC3_JP(X,i,alpha,numFig)

% disp('..::MonteCarlo::..')
[F C] = size(X);
% fprintf('Cells: %d\nFrames: %d\n',C,F);
    
% Observed data
Co = sum(X,2);              % Coactive cells
HO = histc(Co,0:1:C-1);     % Histogram of coactive cells
PO = cumsum(HO)/F;                  % Synchrony probability (OBSERVED)
    
% Simulated data
% fprintf('Simulations: %d\n',i);
for j = 1:i                             % Make i simulations
    for k = 1:C                           
        XS(:,k) = X(randperm(F),k);     % Shuffle spikes for every cell
    end
    CoS = sum(XS,2);                    % Coactive cells (simulated)
    HS = histc(CoS,0:1:C-1);            % Histogram of coactive cells (simulated)
    PS(j,:) = cumsum(HS)/F;             % Synchrony probability (simulated)
end
PE=mean(PS);                            % Mean of Probability (EXPECTED)

% Get propability
P=zeros(1,C);
for j = 1:i
    P = P+double(PO'>=PS(j,:));
end
P = 1-P/i;

% Set threshold from Probability distribution < alpha
th=find(P<alpha);
interval=[0:length(P)-1];

% Probability of get Observed data
figure(numFig)
plot(interval,PO(interval+1),'-r')
title('Probability of obtaining observed data')
hold on
plot(interval,PE(interval+1),'-b')
plot(interval,PO(interval+1)-PE(interval+1)','-g')
plot(interval,P(interval+1)','-M')
plot(interval(th),P(th)+.05,'*m')
hold off
legend('Observed', strcat('Expected (',num2str(i),' simulations)'),...
    'Substract ( Observed - Expected )','Probability (Observed > Expected)')

