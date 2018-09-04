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
%     %% Plot Results 1/3 *************************************
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
        
%         %% Plot Results 2/3 *************************************
%         % WAVES POINTS;
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
                    disp('>>Detrending Linearly ...')
                    mslope=(xxtr(end)-xxtr(1))/length(xxtr);
                    xlinc=mslope*([1:length(xxtr)]-1)+xxtr(1);
                    xdets=xxtr-xlinc;
                    if numel(xdets)>3
                        [~,setPoint]=findpeaks(-xdets,'Npeaks',1,'SortStr','descend');
                        if ~isempty(setPoint)
                            if setPoint>numel(xdets)/2
                                mslope=(xdets(setPoint)-xdets(1))/numel(xdets);
                            else
                                mslope=(xdets(end)-xdets(setPoint))/numel(xdets);
                            end
                            xlincsec=mslope*([1:length(xxtr)]-1)+xdets(1);
                            xlinc=xlinc+xlincsec;
                            xdets=xxtr-xlinc;
                        end
                        if numel(xdets(xdets>=0))<numel(xdets(xdets<0))
                            disp('>>Detrending by RLOESS...')
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
                
                %% Check Peak'sAMplitude between Valleys of Detrended Segment
                if length(xxtr-xlinc)>2
                    [NewLOWwavesFrames,NewZeros] = getwavesframes(xdets,noisex);
                    % Take Care of Big Vallyes:
                    BigVallInd=find(-xdets(NewLOWwavesFrames(2:end-1))>1.5*std(noisex))+1;
                    if ~isempty(BigVallInd)
                        disp('>> Erasing Big Valleys:')
                        NewLOWwavesFrames=NewLOWwavesFrames(setdiff(1:numel(NewLOWwavesFrames),BigVallInd));
                        disp('>>Done.')
                    else
                        disp('- No Negative Distortions -')
                    end
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
            disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Fixing Distortion')
            %% Posterioir Work Arounds
            xdupdate=xd-xlin;
            % xdupdate=xdupdate - var(noisex); % Remove Offset Introduced by Valley-Linear Trending
            [xdenoised,noisex]=mini_denoise(xdupdate);
            % Special Fixes
            % Finale (de)-Bleaching *******************************
            FixS=1;
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
            if max(xdenoised)<std(noisex)
                disp('>> JUST NOISE ')
                xdupdate=xdupdate-xdenoised;
                xdenoised(:)=0; % make it zeros...
            else
                [AmpPEaks,FramPeaks]=findpeaks(xdenoised);
                if isempty(AmpPEaks(AmpPEaks>std(noisex)))
                    disp('>> Fluorescence without Ca++ Transients')
                    xdupdate=xdupdate-xdenoised;
                    xdenoised(:)=0; % make it zeros...
                else
                    if numel(AmpPEaks(AmpPEaks>std(noisex)))==1
                        % Check if Peak Rises FASTER than it Falls
                        FramPeaks=FramPeaks(AmpPEaks>std(noisex));
                        [AmpValls,FramValls]=findpeaks(-xdenoised);
                        if isempty(AmpValls)
                            % No valleys: strange
                            disp('>> Fluorescence without Valleys')
                            xdupdate=xdupdate-xdenoised;
                            xdenoised(:)=0; % make it zeros...
                        else
                            naux=1;
                            while and(and(xdenoised(FramPeaks(1)-naux)>0,~ismember(FramPeaks(1)-naux,FramValls)),FramPeaks(1)-naux>1)
                                naux=naux+1;
                            end
                            RiseN=naux;
                            naux=1;
                            while and(and(xdenoised(FramPeaks(1)+naux)>0,~ismember(FramPeaks(1)+naux,FramValls)),FramPeaks(1)+naux<numel(xdenoised))
                                naux=naux+1
                            end
                            FallN=naux;
                            if RiseN>FallN
                                disp('>> Peak without Ca++ Transients')
                                xdupdate=xdupdate-xdenoised;
                                xdenoised(:)=0; % make it zeros...
                            else
                                disp('>> THERE MIGHT BE Ca++ Transients')
                            end
                        end
                    else
                        disp('>> THERE MIGHT BE Ca++ Transients')
                    end
                    
                end
                
            end
%             % Basal Activity Detectos
%             [~,maxPeakFrame]=findpeaks(xdupdate,'Npeaks',1,'SortStr','descend');
%             if ~isempty(maxPeakFrame)
%                 % Pre Peak
%                 if min(xd(1:maxPeakFrame))>std(noisex)
%                     disp('Possible Initial Basal Activity')
%                 % Post Peak
%                 elseif min(xd(maxPeakFrame:end))>std(noisex)
%                     disp('Possible Final Basal Activity')
%                 end
%                     
%             end
            % check out #######################
            %         plot(xdupdate); pause;
            XDupdate(i,:)=xdupdate;
        end
    end
    %% FEATURE EXTRACTION ############################################
    [SkewSignal(i),SkewNoise(i),SNRbyWT(i),ABratio(i)]=feature_extraction(xdenoised,noisex);
    Xest(i,:)=xdenoised; % SAVE DENOISED  * * * * ** $ $ $ $        OUTPUT
%     %% CHECK RESULTS 3/3
%     h2=subplot(212);
%     plot(XDupdate(i,:)); 
%     hold on;
%     plot([0,numel(XDupdate(i,:))],[0,0],'-.k');
%     plot(xdenoised); 
%     plot([0,numel(xd)],[std(noisex),std(noisex)],'-.r');
%     plot([0,numel(xd)],-[std(noisex),std(noisex)],'-.r');
%     hold off;
%     axis tight; grid on;
%     linkaxes([h1,h2],'x')
%     disp('check')
end % MAIN LOOP    


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