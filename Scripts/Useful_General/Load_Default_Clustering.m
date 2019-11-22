%% Clustering Setups
% SimmEnsTH=0.8;          % Ensemble simmilarity 
MaxEnsembles=10;
TimeDommEns=55;         % Portion of presence of one ensemble
parseq=3/4;             % Portion in a ensemble is activated
RatioHebb=10;           % Minimum % Of the Hebbian Ratio
ActiveNeuronsRatio=0.75;% Ratio of Active Neurons
ActiveTime=0.10;         % Active Time (CAG>0)
% Setup
% alphavalCAG=0.85;       % For CAG  Threshold
SimMethod='hamming';  

% Smoothing CAG for Neural Ensemble Intervals: 
% CONSECUTIVE & INTERLEAVED instances
SmoothWin=fs*30; %            half minute smoothing

% From Features_Condition:
% TypeCycles(1)->Simple
% TypeCycles(2)->Closed
% TypeCycles(3)->Open