% Default Values for Signal Processing ***************************
L=60;              % seconds  of the Fluorophore Response

% AutoRegressive Process Estimated by 'lpc' function *************
p=3;               % AR(p) initial AR order   

% Using fminsearch algorithm to fit AR  **************************
% process impulse response to biexponential Function *************
taus_0= [.75,2,1]; % starting values

% Maximum Lambda Value
MaxLambdaMax=10;

% Theshold for VChS Algorithm
SpuriousTh=1.75;        % Threshold for Peaks