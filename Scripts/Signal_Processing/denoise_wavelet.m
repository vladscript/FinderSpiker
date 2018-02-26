% Funtion to denoise Fluorescence Signals
% By Wavelet Analysis under asumption: noise is iid
% argmin_dbw  { AUC(noise) }
% Input
%   XD:     Detrended Signal
%   X:      Raw Signals
% Output
%   Xest:       Estimated Denoised Signals
%   SNRbyWT:    SNR by wavelet denoising
%   SkewSignal: Skewness of Denoised Signal
%   ABratio:    Max Peak-Valley Amplitude Ratio in Dneoised Signal
%   SkewNoise:  Skewness of Noise
function [Xest,SNRbyWT,SkewSignal,ABratio,SkewNoise,XDupdate]=denoise_wavelet(XD)
%% Setup ****************************
[Ns,Frames]=size(XD);
Xest=zeros(Ns,Frames);
SkewSignal=zeros(Ns,1);
SkewNoise=zeros(Ns,1);
ABratio=zeros(Ns,1);
SNRbyWT=zeros(Ns,1);
XDupdate=XD;
%% FOR---------------------------------------------------------------------
for i=1:Ns
    xd=XD(i,:);
    disp(['------------------>    Denoising Signal: ',num2str(i),'/',num2str(Ns)]);
    % SNR by Wavelet Processing and Autocorrelation Diminishment
    [xdenoised,noisex]=mini_denoise(xd);
    %% DETRENDING FIXING #############################################
    [ValleAMP,ValleN]=findpeaks(-xdenoised);    % Get Valleys of Denoised Signal
    ValleAMPabs=abs(ValleAMP);
    ValleAMPflip=-ValleAMP;
    [pV,binV]=ksdensity(ValleAMPabs);              % pdf of Valley Amplitudes
    [Vp,Vbin,Vwidth]=findpeaks(pV,binV);
    if numel(Vp)>1
        % Take only small amplitudes
        [~,indxsmallAMps]=min(Vbin);
        LOWwaves=ValleAMPflip(ValleAMPflip<=Vbin(indxsmallAMps)+Vwidth(indxsmallAMps));
        LOWwavesFrames=ValleN(ValleAMPflip<=Vbin(indxsmallAMps)+Vwidth(indxsmallAMps));
    else
        LOWwaves=ValleAMP;
        LOWwavesFrames=ValleN;
    end
    LOWwavesFrames=[1,LOWwavesFrames,Frames+1];
    LOWwavesFrames=unique(LOWwavesFrames);
    xlin=[];
    if numel(LOWwavesFrames)>1
        for n=1:numel(LOWwavesFrames)-1
            xxtr=xdenoised(LOWwavesFrames(n):LOWwavesFrames(n+1)-1);
            mslope=(xxtr(end)-xxtr(1))/length(xxtr);
            xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
            xlin=[xlin,xlinc];
        end
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Fixing Distortion')
        xdupdate=xd-xlin;
        [xdenoised,noisex]=mini_denoise(xd);
        % check out #######################
        %         plot(xdupdate); pause;
        XD(i,:)=xdupdate;
    end
    
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

%% NESTED functions
    function [xdenoised,noisex]=mini_denoise(xd)
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
    end    
end