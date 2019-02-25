% Raster Montecarlo 4 (cummulative & fine)
% Get the probability distribution of cells sincronized from raster
% by Montecarlo method using interval inter-activity.
%
% [th P] = MC2_JP(X, i, alpha, numFig)
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
% PO = synchrony probability (OBSERVED)
% PS = synchrony probability (SIMULATED)
%
% ..:: by Jesús E. Pérez-Ortega ::.. Oct-2011

function [th P PO PS] = MC4_JP(X,i,alpha,numFig)

% disp('..::MonteCarlo::..')
[F C] = size(X);
% fprintf('Cells: %d\nFrames: %d\n',C,F);
    
% Observed data
Co = sum(X,2);              % Coactive cells
HO = histc(Co,0:1:C);     % Histogram of coactive cells
PO = cumsum(HO)/F;                  % Synchrony probability (OBSERVED)
    
% Simulated data

for j = 1:i                             % Make i simulations
    XS=zeros(F,C);
    for k = 1:C
        idx=find(X(:,k));
        isi=diff(idx);
        n_isi=length(isi);
        isis=isi(randperm(n_isi));
        idxs=Integrate_JP([randi(F);isis])';
        for l=1:n_isi+1
            if idxs(l)>F
                idxs(l)=idxs(l)-F;
            end
        end
        XS(idxs,k)=1;     % Shuffle spikes for every cell
    end
    CoS = sum(XS,2);                    % Coactive cells (simulated)
    HS = histc(CoS,0:1:C);            % Histogram of coactive cells (simulated)
    PS(j,:) = cumsum(HS)/F;                     % Synchrony probability (simulated)
end
PE=mean(PS);                            % Mean of Probability (EXPECTED)

% Get propability
for j=1:C+1
    P(j)=length(find(PS(:,j)<=PO(j)));
end
P = 1-P/i;

% Set threshold from Probability distribution < alpha
th=find(P(4:C)<alpha,1,'first')+2;
interval=[0:C];
% Probability of get Observed data
figure(numFig)
plot(interval,PO(interval+1),'o-r')
title('Probability of obtaining observed data')
hold on
plot(interval,PE(interval+1),'o-b')
plot(interval,PO(interval+1)-PE(interval+1)','o-g')
plot(interval,P(interval+1)','o-M')
hold off
legend('Observed', strcat('Expected (',num2str(i),' simulations)'),...
    'Substract ( Observed - Expected )','Probability')
