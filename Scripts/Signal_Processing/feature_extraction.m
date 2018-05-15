% Feature Extraction
% Input
%   xdenoised:    clean signal (by wavelet analysis)
%   noisex:       noise(by wavelet analysis)
% Output
%   SkewSignal: Skewness from denoised signal 
%   SkewNoise:  Skewness from noise 
%   SNRbyWT:    Signal Noise Ratio [dB]
%   ABratio:    Max Peaks Aplitudes Ratio

function [SkewSignal,SkewNoise,SNRbyWT,ABratio]=feature_extraction(xdenoised,noisex)
%% FEATURE EXTRACTION ############################################
Aden=findpeaks(xdenoised,'NPeaks',1,'SortStr','descend');
if isempty(Aden)
    Aden=0;
    disp('..................No (+)Peaks');
end
%     Anoise=var(noisex);
%     
Bden=findpeaks(-xdenoised,'NPeaks',1,'SortStr','descend');
if isempty(Bden)
    Bden=0;
    disp('.................No  (-)Peaks');
end    
%     Aratio(i)=Aden/Anoise;  % PEAKS SIGNAL/NOISE RATIO
ABratio=Aden/Bden;   % max/min RATIO
SkewSignal=skewness(xdenoised,0); % $ $ $  OUTPUT
SkewNoise=skewness(noisex,0); % $ $ $  OUTPUT
% SkewRation sensitive to misdetrending distortion
SNRbyWT=10*log(var(xdenoised)/var(noisex)); % 