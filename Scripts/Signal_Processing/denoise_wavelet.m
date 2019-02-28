% Funtion to denoise Fluorescence Signals
% By Wavelet Analysis under asumption: noise is iid
% argmin_dbw  { AUC(noise) }
% & VChS-A: Valley-Checking Segmented Algorithm
% Input
%   XD:         Raw Signal
% Output
%   Xest:       Estimated Denoised Signals
%   SNRbyWT:    SNR by wavelet denoising
%   SkewSignal: Skewness of Denoised Signal
%   ABratio:    Max Peak-Valley Amplitude Ratio in Dneoised Signal
%   SkewNoise:  Skewness of Noise
function [Xest,SNRbyWT,SkewSignal,ABratio,SkewNoise,XDupdate]=denoise_wavelet(XD)
%% Setup ****************************
SpuriousTh=1.25;    % Threshold for Peaks
[Ns,Frames]=size(XD);
Xest=zeros(Ns,Frames);
SkewSignal=zeros(Ns,1);
SkewNoise=zeros(Ns,1);
ABratio=zeros(Ns,1);
SNRbyWT=zeros(Ns,1);
% XDupdate=XD;
%% FOR---------------------------------------------------------------------
for i=1:Ns
    %% GET DATA
    xd=XD(i,:);
    disp(['------------------>    Denoising Signal: ',num2str(i),'/',num2str(Ns)]);
    % SNR by Wavelet Processing and Autocorrelation Diminishment
    [xdenoised,noisex]=mini_denoise(xd);
    %% de-OFFSET
    [~,Valleys]=findpeaks(-xdenoised,'NPeaks',1,'SortStr','descend');
    if isempty(Valleys)
        offset_xd=min(xdenoised)+std(noisex);
    else
        offset_xd=xdenoised(Valleys)+std(noisex);
    end
        
% %     [pxd,binxd]=ksdensity(xdenoised, linspace(min(xdenoised),max(xdenoised),100));
% %     [Ap,Vbin]=findpeaks(pxd,binxd); % modes of pdf valleys
% %     % Caring of too much offset
% %     if numel(Vbin)>1
% %         [~,SortaPeakIndx]=sort(Ap,'descend');
% %         % If more negative mode is of the 2 biigest modes: de -offset
% %         if or(SortaPeakIndx(1)==1,SortaPeakIndx(2)==1)
% %             offset_xd=Vbin(1);
% %         else
% %             offset_xd=Vbin(SortaPeakIndx(end));
% %         end
% %     else
% %         SortaPeakIndx=1;
% %         offset_xd=Vbin;
% %     end
    % if abs(min(xdenoised)-offset_xd)<std(noisex)
    %    offset_xd=Vbin(SortaPeakIndx(1));
    % end
    xdenoised=xdenoised-offset_xd;
    xd=xd-offset_xd;
    XDupdate(i,:)=xd;
    %% Plot Results 1/3 *************************************
%     h1=subplot(211);
%     plot(xd); hold on
%     plot([0,numel(xd)],[std(noisex),std(noisex)],'-.r');
%     plot([0,numel(xd)],-[std(noisex),std(noisex)],'-.r');
%     plot(xdenoised); hold off;
%     axis tight; grid on;
%     disp(i)
    %% CHECK if denoised signal is contained in Standard Deviation of Noise
    if max(xdenoised)<std(noisex)
        % disp('>>JUST NOISE')
        xd=xd-xdenoised;
        xdenoised(:)=0; % make it zeros...
        XDupdate(i,:)=xd;
        % xdenoised(:)=zeros(size())
        % may be ain't happenning ever, man!
    else
        % disp('>> THERE MIGHT BE TRANSIENTS')
        %% DETRENDING FIXING #############################################
        % Get Peaks below noise and waves above noise 
        % Frames would be inflection points to detrend
        [LOWwavesFrames,ZeroCrosses]=getwavesframes(xdenoised,noisex);
        % preLOWwavesFrames=LOWwavesFrames;
        % IGNORE VALLEYS WITH HIGH PROMINENCE!!!
        [~,Framy,~,Promy]=findpeaks(-xdenoised);
        % BigPromyIndx=find(Promy>std(noisex));
        BigValleys=intersect(Framy(Promy>std(noisex)),LOWwavesFrames);
        for nv=1:numel(BigValleys)
            if ismember(BigValleys(nv),LOWwavesFrames)
                Ndel=find(LOWwavesFrames==BigValleys(nv));
                if Ndel>2 && Ndel+1<=numel(LOWwavesFrames)
                    NvecDel=[Ndel-1,Ndel,Ndel+1];
                elseif Ndel-1<2
                    NvecDel=[2,Ndel,Ndel+1];
                elseif Ndel+1>numel(LOWwavesFrames)
                    NvecDel=[1,Ndel,numel(LOWwavesFrames)];
                end
                if NvecDel(end)==numel(LOWwavesFrames)
                    NvecDel=NvecDel(1:2);
                end
                LOWwavesFrames=LOWwavesFrames(setdiff(1:numel(LOWwavesFrames),NvecDel));
            end
        end
        if numel(LOWwavesFrames)==2
            disp('>>Fix of Waves missDetection>>')
            LOWwavesFrames=[LOWwavesFrames(1),numel(xd),LOWwavesFrames(end)];
        elseif numel(LOWwavesFrames)==1
            LOWwavesFrames=[LOWwavesFrames(1),numel(xd),LOWwavesFrames(end)];
        end
            
        % LOWwavesFrames=setdiff(LOWwavesFrames,Framy(Promy>std(noisex)));
        
        %% Plot Results 2/3 *************************************
        % WAVES POINTS;
%         hold on;
%         plot(LOWwavesFrames(2:end-1),xdenoised((LOWwavesFrames(2:end-1))),'*k');
%         hold off
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
                    if LOWwavesFrames(n+2)-1-LOWwavesFrames(n+1)>3
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
                if numel(xxtr)<numel(xd)/3
                    disp('>> Detrending Linearly ...')
                    mslope=(xxtr(end)-xxtr(1))/length(xxtr);
                    xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                    xdets=xxtr-xlinc;
                    if numel(xdets)>3
                        % Biggest Valley:
                        [~,setPoint]=findpeaks(-xdets,'Npeaks',1,'SortStr','descend');
                        if ~isempty(setPoint)
                            % If the Valley is beyond half signal's length
                            if setPoint>numel(xdets)/2
                                % mslope=(xdets(setPoint)-xdets(1))/numel(xdets);
                                mslope=(xxtr(setPoint)-xxtr(1))/numel(xxtr);
                            else
                                mslope=(xdets(end)-xdets(setPoint))/numel(xdets);
                            end
                            % xlincsec=mslope*([1:length(xxtr)]-1)+xdets(1);
                            xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                            % xlinc=xlincsec;
                            xdets=xxtr-xlinc;
                        end
                        if numel(xdets(xdets>=0))<numel(xdets(xdets<0))
                            disp('>> Detrending by RLOESS...')
                            xlinc=smooth(xxtr,numel(xxtr),'rloess')';
                            %xdets=xxtr-xlinc;
                        end
                    end
                else
                    disp('>>Detrending by RLOESS...')
                    xlinc=smooth(xxtr,numel(xxtr),'rloess')';
                end
                xdets=xxtr-xlinc;
                disp('>>...OK')
                
                % CHECK iF it symetric pdf around a mode and ***********
                % and many distributed valleys below zero and noise level
                % Absolute values of Valley are similar o bigger than Peaks
                if numel(xdets)>3
                    Peaks=findpeaks(xdets);
                    Valles=findpeaks(-xdets);
                    if and(~isempty(Peaks),~isempty(Valles))
                        AreAmpsDiff=ttest2(abs(Valles),abs(Peaks));
                        if isnan(AreAmpsDiff)
                            AreAmpsDiff=1;
                        end
                        if ~AreAmpsDiff
                            % xdets=xdets-min(-Valles);
                            xdets=xdets+mean(abs(Valles));
                        end
                    end
                end
                % ******************************************************
                
                %% Check Peak'sAMplitude between Valleys of Detrended Segment
                if numel(xdets)>2
                    [NewLOWwavesFrames,NewZeros] = getwavesframes(xdets,noisex);
                    % Take Care of Big Vallyes NO NEED THANKS 2 GCaMP :
                    % BigVallInd=find(-xdets(NewLOWwavesFrames(2:end-1))>1.5*std(noisex))+1;
                    % if ~isempty(BigVallInd)
                    %     disp('>> Erasing Big Valleys:')
                    %     NewLOWwavesFrames=NewLOWwavesFrames(setdiff(1:numel(NewLOWwavesFrames),BigVallInd));
                    %     disp('>> Done.')
                    % else
                    %     disp('- No Negative Distortions -')
                    % end
                    % NewLOWwavesFrames=NewLOWwavesFrames([1,find(diff(NewLOWwavesFrames(1:end-1))>2)+1,numel(NewLOWwavesFrames)]);
                    if numel(NewLOWwavesFrames)>2
                        xlincnew=[];
                        nn=1;
                        while nn<numel(NewLOWwavesFrames)
                            if nn<=numel(NewLOWwavesFrames)-2 && ...
                                ismember(NewLOWwavesFrames(nn+1),NewZeros)
                                % IF after the Zero-Cross is there a peak
                                if NewLOWwavesFrames(nn+2)-1-NewLOWwavesFrames(nn+1)>3
                                    SecPeak=findpeaks(xdets(NewLOWwavesFrames(nn+1):NewLOWwavesFrames(nn+2)-1));
                                    % Jump Only if there is a peak
                                    if ~isempty(SecPeak)
                                        % n=n+1; % JUMP Zero Cross
                                        NewLOWwavesFrames=[NewLOWwavesFrames(1:nn),NewLOWwavesFrames(nn+2:end)];
                                    % TOO SMALL
                                    elseif NewLOWwavesFrames(nn+1)-NewLOWwavesFrames(nn)<4
                                        NewLOWwavesFrames=[NewLOWwavesFrames(1:nn),NewLOWwavesFrames(nn+2:end)];
                                    end
                                else
                                    % n=n+1; % JUMP Zero Cross
                                    NewLOWwavesFrames=[NewLOWwavesFrames(1:nn),NewLOWwavesFrames(nn+2:end)];
                                end
                            end
                            xlincs=getlinearsegment(xdets(NewLOWwavesFrames(nn):NewLOWwavesFrames(nn+1)-1),std(noisex),nn);
                            xlincnew=[xlincnew,xlincs];
                            % Should I update the 1st samples of the next 
                            % interval as a detrended signal ?????????
                            disp(nn)
                            nn=nn+1;
                        end
                        xlincnew(xlincnew>xdets)=xdets(xlincnew>xdets);
                        xxtr=xdets;
                        xlinc=xlincnew+xlinc;                
                        % AmpPeak=findpeaks(xxtr);   % Peaks
                    else
                        % AmpPeak=findpeaks(xxtr-xlinc);   % Peaks
                        
                    end

                end
                xlin=[xlin,xlinc];
                n=n+1; % NEXT VALLEY-ZeroCross
            end
            disp('~   ~  ~ ~ Fixing Distortion')
            %% Posterioir Work Arounds
            xdupdate=xd-xlin;            
            % xdupdate=xdupdate - var(noisex); % Remove Offset Introduced by Valley-Linear Trending
            [xdenoised,noisex]=mini_denoise(xdupdate);
            
            % Make it Non-Negative: OFFSETTIN SOME SIGNALS
            %xdupdate=xdupdate-(min(xdenoised));
            %xdenoised=xdenoised-(min(xdenoised));
%             % TOO MUCH NEGATIVITY FIXING:
%             if min(xdenoised)<-std(noisex)
%                 deltaoffset=abs(-std(noisex)-min(xdenoised));
%                 xdupdate=xdupdate+deltaoffset;
%                 xdenoised=xdenoised+deltaoffset;
%                 disp('>>Fixed Negativity: :D')
%             end
            
            % Final Detrending
            [LOWwavesFrames,~]=getwavesframes(xdenoised,noisex);
            NEGpoints=find(xdenoised(LOWwavesFrames(2:end-1))<0);
            if ~isempty(NEGpoints)
                delineframes=LOWwavesFrames(NEGpoints+1);
                disp('>> Detrending: [final fix]')
                xdeline=getdeline(delineframes,xdenoised);
                xdupdate=xdupdate-xdeline;
                xdenoised=xdenoised-xdeline;
            end
            disp('>> Detrending: [OK]')
            
            % Special Fixes
            % Finale (de)-Bleaching *******************************
            % FixS=1;
            if xdenoised(end)==max(xdenoised) && numel(LOWwavesFrames)>3
                dxden=diff(xdenoised);
                FixS=LOWwavesFrames(end-1);
                xdupdate(end-FixS:end)=xdupdate(end-FixS:end)-xdenoised(end-FixS:end);
                xdenoised(end-FixS:end)=0;
                disp('>>>> Bleaching Fixed')
            end
            % Check Inter Points
            NP=LOWwavesFrames(2:end-1);
            disp('Check Point .')
            for k=1:numel(NP)
                fprintf('.');
                preAmp=xdenoised(NP(k));
                if NP(k)>=numel(xdenoised)
                    posAmp=xdenoised(NP(k));
                else
                    posAmp=xdenoised(NP(k)+1);
                end
                if abs(preAmp-posAmp)>std(noisex)
                    disp('ups...')
                end
            end
            fprintf('\n');
            if max(xdenoised)<SpuriousTh*std(noisex)
                disp('>> JUST NOISE ')
                xdupdate=xdupdate-xdenoised;
                xdenoised(:)=0; % make it zeros...
            else
                [AmpPEaks,FramPeaks]=findpeaks(xdenoised);
                FramPeaks=FramPeaks(AmpPEaks>SpuriousTh*std(noisex));
                [AmpValls,FramValls]=findpeaks(-xdenoised);
                nAmps=find(AmpPEaks>SpuriousTh*std(noisex));
                if isempty(AmpPEaks(AmpPEaks>SpuriousTh*std(noisex)))
                    disp('>> Fluorescence without Ca++ Transients')
                    xdupdate=xdupdate-xdenoised;
                    xdenoised(:)=0; % make it zeros...
                else
                    % For the big&single slow yummies *************************
                    if numel(AmpPEaks(AmpPEaks>SpuriousTh*std(noisex)))==1
                        if isempty(AmpValls)
                            % No valleys: strange
                            disp('>> Fluorescence without Valleys')
                            xdupdate=xdupdate-xdenoised;
                            xdenoised(:)=0; % make it zeros...
                        else
                            fn=1;
                            naux=1;
                            while and(and(xdenoised(FramPeaks(fn)-naux)>0,~ismember(FramPeaks(fn)-naux,FramValls)),FramPeaks(fn)-naux>1)
                                naux=naux+1;
                            end
                            RiseN=naux;
                            naux=1;
                            while and(and(xdenoised(FramPeaks(fn)+naux)>0,~ismember(FramPeaks(fn)+naux,FramValls)),FramPeaks(fn)+naux<numel(xdenoised))
                                naux=naux+1;
                            end
                            FallN=naux;
                            if RiseN>FallN
                                % xdupdate(FramPeaks-RiseN:FramPeaks-FallN)=xdupdate(FramPeaks-RiseN:FramPeaks-FallN)-xdenoised(FramPeaks-RiseN:FramPeaks-FallN);
                                % xdenoised(FramPeaks-RiseN:FramPeaks-FallN)=0;
                                xdupdate=xdupdate-xdenoised;
                                xdenoised(:)=0; % make it zeros...
                                disp('>> Small Peak without Ca++ Transients')
                            else
                                disp('>> THERE MIGHT BE Ca++ Transients')
                            end
                        end
                    else
                        OkPeaks=zeros(numel(FramPeaks),1);
                        for fn=1:numel(FramPeaks)
                            naux=1;
                            while and(and(xdenoised(FramPeaks(fn)-naux)>0,~ismember(FramPeaks(fn)-naux,FramValls)),FramPeaks(fn)-naux>1)
                                naux=naux+1;
                            end
                            RiseN=naux;
                            naux=1;
                            while and(and(xdenoised(FramPeaks(fn)+naux)>0,~ismember(FramPeaks(fn)+naux,FramValls)),FramPeaks(fn)+naux<numel(xdenoised))
                                naux=naux+1;
                            end
                            FallN=naux;
                            if FallN<RiseN
                                disp('>> Symetric Waveform.')
                                % xdupdate=xdupdate-xdenoised;
                                % xdenoised(:)=0; % make it zeros...
                            else
                                disp('>> THERE MIGHT BE Ca++ Transients')
                                OkPeaks(fn)=1;
                            end
                        end
                        if sum(OkPeaks==0)==numel(FramPeaks)
                            xdupdate=xdupdate-xdenoised;
                            xdenoised(:)=0; % make it zeros...
                            disp('>> All the Small PEAKS without Ca++ Transients')
                        end
                    end
                    
                end
                
            end
            % Update Value:
            XDupdate(i,:)=xdupdate;
        end
    end
    %% FEATURE EXTRACTION ############################################
    [SkewSignal(i),SkewNoise(i),SNRbyWT(i),ABratio(i)]=feature_extraction(xdenoised,noisex);
    Xest(i,:)=xdenoised; % SAVE DENOISED  * * * * ** $ $ $ $        OUTPUT
 %%   CHECK RESULTS 3/3
%     h2=subplot(212);
%     plot(XDupdate(i,:)); 
%     hold on;
%     plot(xdenoised,'LineWidth',2); 
%     plot([0,numel(XDupdate(i,:))],[0,0],'-.k');
%     plot([0,numel(xd)],[std(noisex),std(noisex)],'-.g');
%     plot([0,numel(xd)],-[std(noisex),std(noisex)],'-.g');
%     hold off;
%     axis tight; grid on;
%     title(num2str(i))
%     linkaxes([h1,h2],'x')
%     disp('check')
end % MAIN LOOP    %     




%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
end % END OF THE WORLD####################################################


function [samplestodetrend,ylinv1]=lineardetbigmode(binx,ModesX,x)
    % Get Samples belonging from least sample to the biggest mode
    Frames=length(x);
    ModePos=find(binx==ModesX(1));
    % dpx=diff(px);
    samplestodetrend= intersect( find(x>min(x)),find(x<binx(ModePos)) );
    [cl,~,mucl] = polyfit(samplestodetrend,x(samplestodetrend),1);
    ylinv1=polyval(cl,1:Frames,[],mucl); % Linear Trending
end


function yss = smoothbysegments(x,y,NC)
    [~,PeaksY]=findpeaks(y);
    [~,ValleY]=findpeaks(-y);
    Curves=[makerowvector(PeaksY),makerowvector(ValleY)];
    if numel(Curves)>NC
        disp('   -^-v- Curvy Trending )))')
        Curves=sort(Curves);
        frameA=1;
        yss=zeros(size(x));
        for n=1:numel(Curves)+1
            if n<=numel(Curves)
                frameB=Curves(n);
            else
                frameB=numel(x);
            end
            xseg=x(frameA:frameB);
            yseg=smooth(xseg,numel(xseg),'rloess'); 
            % (RE)-CURSIVENESS:############ NOT YET
            % yseg=smoothbysegments(xseg,yseg);
            yss(frameA:frameB)=yseg;
            frameA=frameB;
        end
        
    else
        yss=y';
    end
end