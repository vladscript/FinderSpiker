% Funtion to denoise Fluorescence Signals
% By Wavelet Analysis under asumption: noise is iid
% argmin_dbw  { AUC(noise) }
% Input
%   XD:     Detrended Signal
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
    %% de-OFFSET
    [pxd,binxd]=ksdensity(xdenoised, linspace(min(xdenoised),max(xdenoised),100));
    [Ap,Vbin]=findpeaks(pxd,binxd); % modes of pdf valleys
    % Caring of too much offset
    if numel(Vbin)>1
        [~,SortaPeakIndx]=sort(Ap,'descend');
        % If more negative mode is of the 2 biigest modes: de -offset
        if or(SortaPeakIndx(1)==1,SortaPeakIndx(2)==1)
            offset_xd=Vbin(1);
        else
            offset_xd=Vbin(SortaPeakIndx(1));
        end
    else
        SortaPeakIndx=1;
        offset_xd=Vbin;
    end
    if abs(min(xdenoised)-offset_xd)<std(noisex)
        offset_xd=Vbin(SortaPeakIndx(1));
    end
    xdenoised=xdenoised-offset_xd;
    xd=xd-offset_xd;
    %% DETRENDING FIXING #############################################
    [ValleAMP,ValleN]=findpeaks(-xdenoised);    % Get Valleys of Denoised Signal
    if ~isempty(ValleAMP)
        ValleAMPabs=abs(ValleAMP);
        ValleAMPflip=-ValleAMP;              % actual values
        [pV,binV]=ksdensity(xdenoised(ValleN));    % pdf of Valley Amplitudes
        [Vp,Vbin,Vwidth]=findpeaks(pV,binV); % modes of pdf valleys
        if numel(Vp)>0
            % Take only small amplitudes
            [~,indxsmallAMps]=max(Vp);
            ThAmpValley=max(Vbin(indxsmallAMps)+Vwidth(indxsmallAMps),std(noisex));
            LOWwaves=ValleAMPflip(ValleAMPflip<=ThAmpValley);
            LOWwavesFrames=ValleN(ValleAMPflip<=ThAmpValley);
        else
            %if ~isempty(ValleAMPabs)
            % ksdensity error-> use histogram
            [Counts,Edges]=histcounts(xdenoised(ValleN),length(ValleAMP));
            [~,indxbin]=max(Counts);
            LOWwaves=ValleAMPflip(ValleAMPflip<=Edges(indxbin+1));
            LOWwavesFrames=ValleN(ValleAMPflip<=Edges(indxbin+1));
            disp('by histo[...]gram')
            %else
            %    LOWwaves=ValleAMP;
            %    LOWwavesFrames=ValleN;
            %end
        end
    else
        LOWwaves=[];
        LOWwavesFrames=[];
        disp('-------------------------No distortion issues')
    end
    
    LOWwavesFramesAll=unique(LOWwavesFrames);
    LOWwavesFrames=LOWwavesFramesAll(LOWwaves<std(noisex));
    LOWwavesFrames=[1,LOWwavesFrames,Frames+1]; % ok that '+1'
    ZeroCrosses=find(diff(sign(xdenoised)));    % Zero Crosses
    LOWwavesFrames=unique([LOWwavesFrames,ZeroCrosses]);
    % Ignore low spaced waves
    LOWwavesFrames=LOWwavesFrames(setdiff (1:length(LOWwavesFrames), find(diff(LOWwavesFrames)<2)+1) );
    ZerosinWavesIndx=find(ismember(LOWwavesFrames,ZeroCrosses));
    IgnoreWaves=find(diff(ZerosinWavesIndx)==1)+1;
    AcceptWaves=setdiff(1:numel(LOWwavesFrames),ZerosinWavesIndx(IgnoreWaves));
    LOWwavesFrames=LOWwavesFrames(AcceptWaves);
%     %% Plot Results *************************************
%     h1=subplot(211);
%     plot(xd); hold on
%     plot(xdenoised);
%     % plot(xlin);
%     plot(LOWwavesFrames(2:end-1),xdenoised((LOWwavesFrames(2:end-1))),'*k');
%     plot([0,numel(xd)],[std(noisex),std(noisex)],'-.r');
%     plot([0,numel(xd)],-[std(noisex),std(noisex)],'-.r');
%     hold off
%     axis tight; grid on;
%     disp(i)
    
    %% AVChS: Algorithm Valley Checking Section *************************************
    xlin=[];
    n=1;
    if numel(LOWwavesFrames)>2
        % LOOP VALLEYs
        while n<numel(LOWwavesFrames)
            % IF there is a follwoing Zero-Cross:
            if n<=numel(LOWwavesFrames)-2 && ...
                    ismember(LOWwavesFrames(n+1),ZeroCrosses)
                % IF after the Zero-Cross is there a peak
                if LOWwavesFrames(n+2)-1-LOWwavesFrames(n+1)>2
                    SecPeak=findpeaks(xdenoised(LOWwavesFrames(n+1):LOWwavesFrames(n+2)-1));
                    if ~isempty(SecPeak)
                        % n=n+1; % JUMP Zero Cross
                        LOWwavesFrames=[LOWwavesFrames(1:n),LOWwavesFrames(n+2:end)];
                    end
                else
                    % n=n+1; % JUMP Zero Cross
                    LOWwavesFrames=[LOWwavesFrames(1:n),LOWwavesFrames(n+2:end)];
                end
            end
            xxtr=xdenoised(LOWwavesFrames(n):LOWwavesFrames(n+1)-1);
            
            % Avoid detrend higher points
            if  and( xxtr(1)>std(noisex) , xxtr(end)>-std(noisex))
                % 1st and Final Sample Test above noise
                if n==1
                    if length(xxtr)>2
                        mslope=(xxtr(end)-xxtr(1))/length(xxtr);
                        xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                        Apeaks=findpeaks(xxtr-xlinc);
                    else
                        Apeaks=[];
                    end
                    if isempty(Apeaks)
                        xexp=fit([1:length(xxtr)]',xxtr','exp1');
                        xlinc=xexp(1:length(xxtr))';
                    else
                        mslope=0;
                        xlinc=mslope*([1:length(xxtr)]-1);
                    end
                    % Pure Decaying at the Start
                else
                    mslope=0;
                    xlinc=mslope*([1:length(xxtr)]-1);
                    % Possible Calcium Transient
                end
            elseif and(xxtr(1)<-std(noisex),xxtr(end)<-std(noisex))
                % 1st Sample below Noise & Final Sample below noise
                %if and( numel(xxtr(xxtr>max([xxtr(1),xxtr(end)])))>numel(xxtr(xxtr<max([xxtr(1),xxtr(end)]))),...
                %    or(xxtr(1)>std(noisex),xxtr(end)>-std(noisex)) )
                if numel(xxtr(xxtr>max([xxtr(1),xxtr(end)])))>numel(xxtr(xxtr<max([xxtr(1),xxtr(end)])))
                    mslope=(xxtr(end)-xxtr(1))/length(xxtr);
                    if  mslope>0
                        mslope=0;
                    end
                    xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                    disp('Ca2+ Transient')
                else
                    xlinc=xxtr;
                end
            elseif and(xxtr(1)<-std(noisex),xxtr(end)>-std(noisex))
                % 1st Sample below & Final Sample above noise
                if or( numel(xxtr(xxtr>std(noisex)))>numel(xxtr(xxtr<0)),...
                    or(xxtr(1)>std(noisex),xxtr(end)>-std(noisex)) )
                    mslope=(xxtr(end)-xxtr(1))/length(xxtr);
                    xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                    disp('Ca2+ Transient')
                else
                    xlinc=xxtr;
                end
            else
                if xxtr(end)>std(noisex)
                    mslope=0;
                    xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                    disp('Ca2+ Transient')
                else
                    % Otherwise Line
                    mslope=(xxtr(end)-xxtr(1))/length(xxtr);
                    xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                    if sum(xlinc>xxtr)>numel(xxtr)/2
                        % If line is above data, use smooth:
                        xlinc=xxtr;
                    end
                end
            end
            
            % Check Peak'sAMplitude between Valleys of Detrended Segment
            if length(xxtr-xlinc)>2
                AmpPeak=findpeaks((xxtr-xlinc),'SortStr','descend');    % Peaks
                % AmpValley=findpeaks(-(xxtr-xlinc),'SortStr','descend'); % Valleys
                if and(~isempty(AmpPeak),and(n>1,n<numel(LOWwavesFrames)-1))
                    if AmpPeak(1)>=std(noisex) % PEAK ABOVE NOISE: linear
                        disp('-OK: Ca Transient-')
                        mslope=(xxtr(end)-xxtr(1))/length(xxtr);
                        xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                    else
                        xlinc=xxtr;            % PEAK BELOW NOISE : smooth
                        disp('low movements <A>')
                    end
                else                           % NO PEAKS: linear
                    disp('low movements <B>')
                end
            end
            xlin=[xlin,xlinc];
            n=n+1; % NEXT VALLEY-ZeroCross
        end
        disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Fixing Distortion')
        xdupdate=xd-xlin;
        % xdupdate=xdupdate - var(noisex); % Remove Offset Introduced by Valley-Linear Trending
        [xdenoised,noisex]=mini_denoise(xdupdate);
        % check out #######################
        %         plot(xdupdate); pause;
        XDupdate(i,:)=xdupdate;
    end
%     % CHECK
%     h2=subplot(212);
%     plot(XDupdate(i,:)); 
%     hold on;
%     plot([0,numel(XDupdate(i,:))],[0,0],'-.r'); hold off;
%     axis tight; grid on;
%     linkaxes([h1,h2],'x')
%     disp('dick')
    %% FEATURE EXTRACTION ############################################
    [SkewSignal(i),SkewNoise(i),SNRbyWT(i),ABratio(i)]=feature_extraction(xdenoised,noisex);
    Xest(i,:)=xdenoised; % SAVE DENOISED  * * * * ** $ $ $ $        OUTPUT
    
    
end

end