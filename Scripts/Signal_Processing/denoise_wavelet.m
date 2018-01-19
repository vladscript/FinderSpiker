% Funtion to denoise Fluorescence Signals
% By Wavelet Analysis under asumption: noise is iid
% argmin_dbw  { AUC(noise) }
% Input
%   XD:     Detrended Signal
% Output
%   Xest:       Estimated Denoised Signals
%   SNRbyWT:    SNR by wavelet denoising
%   SkewSignal: Skewness of Denoised Signal
%   ABratio:    Max Peak-Valley Amplitude Ratio
%   SkewNoise:  Skewness of Noise
function [Xest,SNRbyWT,SkewSignal,ABratio,SkewNoise]=denoise_wavelet(XD)
[Ns,Frames]=size(XD);
Xest=zeros(Ns,Frames);
SkewSignal=zeros(Ns,1);
SkewNoise=zeros(Ns,1);
ABratio=zeros(Ns,1);
SNRbyWT=zeros(Ns,1);
%% FOR---------
for i=1:Ns
    xd=XD(i,:);
    disp(['------------------>    Denoising Signal: ',num2str(i),'/',num2str(Ns)]);
    % SNR by Wavelet Processing and Autocorrelation Diminishment
    dbw=1;          % Degree (Order )Wavelet
    gradacf=-10;
    acf_pre=1;
    ondeleta='sym8'; % Sort of Wavelet
    while ~or(gradacf>=0,dbw>9)
        disp(['DENOISING by wavelet level . . . . . ',num2str(dbw)])
        %         [xdenoised,~,~,~,~] = cmddenoise(xd,ondeleta,dbw);
        %          xdenoised=smooth(xd,SWS,'rloess'); ALTERNATIVE
        xdenoised=waveletdenoise(xd,ondeleta,dbw);
        noise=xd-xdenoised;
        acf=autocorr(noise,1);
        acf_act=abs(acf(end));
        gradacf=acf_act-acf_pre;
        acf_pre=acf_act;    % update ACF value
        dbw=dbw+1;          % Increase wavelet level
    end
    disp(' denoised <<<<<<<<<<<<<<<<<<<<<<<')
    dbw_opt=dbw-2; % Beacuse it detected +1 and increased it anyway +1
    %     [xdenoised,~,~,~,~] = cmddenoise(xd,ondeleta,dbw_opt);
    xdenoised=waveletdenoise(xd,ondeleta,dbw_opt);
    noisex=xd-xdenoised;
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
    ABratio(i)=Aden/Bden;   % max/min RATIO
    SkewSignal(i)=skewness(xdenoised,0); % $ $ $  OUTPUT
%     if isnan(SkewSignal(i))
%         disp('ass')
%     end
    SkewNoise(i)=skewness(noisex,0); % $ $ $  OUTPUT
    % SkewRation sensitive to misdetrending distortion
    
    SNRbyWT(i)=10*log(var(xdenoised)/var(noise)); % $ $ $ $ $       OUTPUT
    Xest(i,:)=xdenoised; % SAVE DENOISED  * * * * ** $ $ $ $        OUTPUT
end