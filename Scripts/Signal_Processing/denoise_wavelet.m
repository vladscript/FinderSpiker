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
    % Get Peaks below noise and waves above noise 
    % Frames would be inflection points to detrend
    [LOWwavesFrames,ZeroCrosses]=getwavesframes(xdenoised,noisex);
%     %% Plot Results 1/2 *************************************
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
    
    %% VChS-A: Valley-Checking Segmented Algorithm *************************************
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
            
            %% DEFINE XLINC (detrendig by parts)
            
            xlinc=getlinearsegment(xxtr,noisex,n);
            %% Check Peak'sAMplitude between Valleys of Detrended Segment
            if length(xxtr-xlinc)>2
                [NewLOWwavesFrames,NewZeros] = getwavesframes(xxtr-xlinc,noisex);
                if numel(NewLOWwavesFrames)>2
                    xlincnew=[];
                    xdets=xxtr-xlinc;
                    nn=1;
                    while nn<numel(NewLOWwavesFrames)
                        if nn<=numel(NewLOWwavesFrames)-2 && ...
                            ismember(NewLOWwavesFrames(nn+1),NewZeros)
                            % IF after the Zero-Cross is there a peak
                            if NewLOWwavesFrames(nn+2)-1-NewLOWwavesFrames(nn+1)>2
                                SecPeak=findpeaks(xdets(NewLOWwavesFrames(nn+1):NewLOWwavesFrames(nn+2)-1));
                                if ~isempty(SecPeak)
                                    % n=n+1; % JUMP Zero Cross
                                    NewLOWwavesFrames=[NewLOWwavesFrames(1:nn),NewLOWwavesFrames(nn+2:end)];
                                end
                            else
                                % n=n+1; % JUMP Zero Cross
                                NewLOWwavesFrames=[NewLOWwavesFrames(1:nn),NewLOWwavesFrames(nn+2:end)];
                            end
                        end
                        xlincs=getlinearsegment(xdets(NewLOWwavesFrames(nn):NewLOWwavesFrames(nn+1)-1),noisex,nn);
                        xlincnew=[xlincnew,xlincs];
                        disp(nn)
                        nn=nn+1;
                    end                
                xxtr=xdets;
                xlinc=xlincnew+xlinc;
                AmpPeak=findpeaks(xxtr);   % Peaks
                else
                    AmpPeak=findpeaks(xxtr-xlinc);   % Peaks
                end
                
                if and(~isempty(AmpPeak),and(n>1,n<numel(LOWwavesFrames)-1))
                    if max(AmpPeak)>=std(noisex) % PEAK ABOVE NOISE: linear
                        disp('-OK: Ca Transient-')
                        % mslope=(xxtr(end)-xxtr(1))/length(xxtr);
                        % xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                    else
                        % xlinc=xxtr;            % PEAK BELOW NOISE : smooth
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
        %% Posterioir Work Arounds
        xdupdate=xd-xlin;
        % xdupdate=xdupdate - var(noisex); % Remove Offset Introduced by Valley-Linear Trending
        [xdenoised,noisex]=mini_denoise(xdupdate);
        % Special Fixes
        % Initial Bleaching *******************************
%         FixS=1;
%         if xdenoised(1)>std(noisex)
%             dxden=diff(xdenoised);
%             while dxden(FixS)<0
%                 FixS=FixS+1;
%             end
%             if xdenoised(FixS)<-std(noisex)
%                while xdenoised(FixS)<0
%                     FixS=FixS+1;
%                end
%             end
%             xdupdate(1:FixS)=xdupdate(1:FixS)-xdenoised(1:FixS);
%             xdenoised(1:FixS)=0;
%             disp('>>>> Bleaching Fixed')
%         end
% Finale (de)-Bleaching *******************************
        FixS=1;
        if xdenoised(end)==max(xdenoised)
            dxden=diff(xdenoised);
            FixS=LOWwavesFrames(end-1);
            xdupdate(end-FixS:end)=xdupdate(end-FixS:end)-xdenoised(end-FixS:end);
            xdenoised(end-FixS:end)=0;
            disp('>>>> Bleaching Fixed')
        end
        % check out #######################
        %         plot(xdupdate); pause;
        XDupdate(i,:)=xdupdate;
    end
% %     % CHECK 2/2
%     h2=subplot(212);
%     plot(XDupdate(i,:)); 
%     hold on;
%     plot([0,numel(XDupdate(i,:))],[0,0],'-.r');
%     plot(xdenoised); 
%     plot([0,numel(xd)],[std(noisex),std(noisex)],'-.r');
%     plot([0,numel(xd)],-[std(noisex),std(noisex)],'-.r');
%     hold off;
%     axis tight; grid on;
%     linkaxes([h1,h2],'x')
%     disp('dick')
    %% FEATURE EXTRACTION ############################################
    [SkewSignal(i),SkewNoise(i),SNRbyWT(i),ABratio(i)]=feature_extraction(xdenoised,noisex);
    Xest(i,:)=xdenoised; % SAVE DENOISED  * * * * ** $ $ $ $        OUTPUT
    
    
end

end