% Funtion that estimates the coefficients of Auto Regressive Process
% Input
%   x: Signal
%   p: Order
%   fs: sampling frequency
%   L: Duration of Response in Seconds
%   taus_0: initial taus
% Output
%   r: Biexponential Response Fit
%   ARc: Coefficients
%   TAUS: Biexponential Paramters
function [FR,ARcoefficients,TAUS]=AR_Estimation(XD,p,fs,L,tau_0)
% Setup
samplesL=round(L*fs);
[C,~]=size(XD); % Cells & Frames
impulse    = zeros(samplesL,1); impulse(1) = 1;    % Impulse
ARcoefficients=[];
TAUS=[];
FR=[];
% Main Loop
for c=1:C
    x=XD(c,:);
    % Estimate Coefficients
    ARc = lpc(x,p);   % coefficients of IIR filter    
    ARcoefficients=[ARcoefficients;ARc];
    r = filter(1,ARc,impulse);    % Response FIR (inverse filter)
    [~,poStart]=findpeaks(-r,'NPeaks',1);
    if isempty(poStart)
        poStart=1;
    end
    Amp_offset=r(poStart);
    r(1:poStart)=0;
    % r(1)=0;
    % Biexponential Parameters
    ts=1/fs;
    t_r = linspace(0,length(r)*ts,length(r))';   
    [taus,~] = lptobiexp(tau_0,t_r,r);  % find taus
    TAUS=[TAUS;taus];
    r_bi = taus(3)*(exp(-t_r/taus(2)) - exp(-t_r/taus(1)));
    r = r_bi; % use bi-exponential vector
    FR=[FR;r'];
end