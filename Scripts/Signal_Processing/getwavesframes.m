% Function to get frames in waves that are above or below
% Standard Divation of Noise
% Input
%   xdenoised: denoised signal (smoothed)
%   noise: noise array (xoriginal-xdenoised)
% Output
%   Frames of Waves: Peaks (above noise) and Valleys (below noise)

function [LOWwavesFrames,ZeroCrosses] = getwavesframes(xdenoised,noisex)
    Frames=numel(xdenoised);
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

end