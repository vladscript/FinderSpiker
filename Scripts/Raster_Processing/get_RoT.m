% Function that Detects Calcium Transients
% from onset spikes and denoised fluorescence signal
% Input
%   R:      Raster Matrix: Cells x Frames
%   XEST:   Denoised Fluo Signals
% Output
%   NoT     Number of Transient per Cell
%   ITI     Inter Transient Interval    units:[samples]
%   LT      Latency Transient           units:[samples]
function [NoT,ITI,LT]=get_RoT(R,XEST)
% Setup
NC=size(R,1);
NoT=zeros(NC,1);
ITI=[];
LT=[];
% Main Loop
for c=1:NC
    % data
    r=R(c,:);       % Activity Row
    xe=XEST(c,:);   % denoised Signal
    if max(xe)~=min(xe)
        ThAmp=0.25*(max(xe)-min(xe));
    else
        ThAmp=0;
    end
    % Detect Peaks And Valleys
    [~,Peaks]=findpeaks(xe,'MinPeakHeight',ThAmp);
    [~,Onsets]=findpeaks(-xe);   
    if ~isempty(Peaks)&&~isempty(Onsets)
        if Peaks(1)<Onsets(1)
            Onsets=[1,Onsets];
        end
        if Onsets(end)>Peaks(end)
            Peaks=[Peaks,numel(xe)];
        end
    end
    
%     % Checking Stuff
%     ax1=subplot(211);plot(xe)
%     hold on;
%     plot(Peaks,xe(Peaks),'k*')
%     plot(Onsets,xe(Onsets),'k^')
%     hold off;
%     axis tight; grid on;
%     ax2=subplot(212);plot(r)
%     axis tight; grid on;
%     linkaxes([ax1,ax2],'x')
%     NoT(c)
    
    % Valley-Peak Loop
    for no=1:numel(Peaks)
        % ActualOnset=Onsets(no);
        % ActualPeak=Peaks(find(Peaks-ActualOnset>0,1));
        
        ActualPeak=Peaks(no);
        OnsetsIndx=find(Onsets-ActualPeak<0);
        if isempty(OnsetsIndx)
            OnsetsIndx=1;
            Onsets=1;
        end
        ActualOnset=Onsets(OnsetsIndx(end));
        
        if isempty(ActualPeak)
            ActualPeak=numel(xe);
        end
        % Detect Transients in between Valley-Peak
        if sum(r(ActualOnset:ActualPeak))>0
            % Transient Counter
            NoT(c)=NoT(c)+1;
            % Inter Transient Interval
            if NoT(c)>1
                ITI=[ITI;ActualOnset-prePeak];
            end
            prePeak=ActualPeak;
            % Transient Duration
            LT=[LT;ActualPeak-ActualOnset+1];
                     
            fprintf('^-')
        end
    end
   fprintf('\n') 
end