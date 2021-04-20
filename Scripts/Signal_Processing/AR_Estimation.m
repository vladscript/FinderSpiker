% Funtion that estimates the coefficients of Auto Regressive Process
% Input
%   XD:     Matrix of signals: Signals x Samples
%   p:      Order of Linear Predictor Coefficients (FIR-filter)
%   fs:     sampling frequency @ Hz
%   L:      Duration of Response in Seconds
%   taus_0: initial taus of the AR(p) process
% Output
%   FR:      Fitted Biexponential Responses to the AR process
%   ARc:     Coefficients of AR(p)
%   TAUS:    Biexponential Paramters: [tau_rise,tau_fall,gain]
function [FR,ARcoefficients,TAUS]=AR_Estimation(XD,p,fs,L,tau_0)
%% Setup
samplesL=round(L*fs);
impulse    = zeros(samplesL,1); impulse(1) = 1;    % Impulse
ARcoefficients=zeros(size(XD,1),p+1);
TAUS=zeros(size(XD,1),3);
FR=zeros(size(XD,1),samplesL);
%% Main Loop
fprintf('>>AR(%i) processes:',p)
for c=1:size(XD,1)
    x=XD(c,:);
    % Estimate Coefficients
    ARc = lpc(x,p);   % coefficients of IIR filter    
    ARcoefficients(c,:)=ARc;
    r = filter(1,ARc,impulse);    % Response FIR (inverse filter)
    [~,poStart]=findpeaks(-r,'NPeaks',1);
    if isempty(poStart)
        poStart=1;
    end
    % Amp_offset=r(poStart);
    r(1:poStart)=0;
    % r(1)=0;
    % Biexponential Parameters
    ts=1/fs;
    t_r = linspace(0,length(r)*ts,length(r))';   
    [taus,~] = lptobiexp(tau_0,t_r,r);  % find taus
    TAUS(c,:)=taus;
    r_bi = taus(3)*(exp(-t_r/taus(2)) - exp(-t_r/taus(1)));
    r = r_bi; % use bi-exponential vector
    FR(c,:)=r';
    fprintf('%i ',c)
end
fprintf('\n')