% Function that Detects Calcium Transients
% from onset spikes and denoised fluorescence signal
% Input
%   R:      Raster Matrix: Cells x Frames
%   XEST:   Denoised Fluo Signals
% Output
%   RoT     Rate of Transient per Cell
function RoT=get_RoT(R,XEST)
% Setup
NC=size(R,1);
RoT=zeros(NC,1);
for c=1:NC
    r=R(c,:);       % Activity Row
    xe=XEST(c,:);   % denoised Signal
    % Detect Peaks And Valleys
    [~,Peaks]=findpeaks(xe);
    [~,Onsets]=findpeaks(-xe);
    for no=1:numel(Onsets)
        ActualOnset=Onsets(no);
        ActualPeak=Peaks(find(Peaks-ActualOnset>0,1));
        if isempty(ActualPeak)
            ActualPeak=numel(xe);
        end
        % Detect Transients in between Valley-Peak
        if sum(r(ActualOnset:ActualPeak))>0
            RoT(c)=RoT(c)+1;
            fprintf('_|\\_')
        end
    end
    fprintf('\n')
end