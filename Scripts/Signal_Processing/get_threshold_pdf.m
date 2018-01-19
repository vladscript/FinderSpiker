% Funtion to get Positive Threshold of a Feature
% where p(signal)>p(noise)
% using indexes of detected signals and noise
function Th_pdf=get_threshold_pdf(Feature,signalindx,noiseindx)
if ~isempty(noiseindx) % Detected Signal and Noise Records
    supportfeat=linspace(min([min(Feature(signalindx)),min(Feature(noiseindx))]),...
        max([max(Feature(signalindx)),max(Feature(noiseindx))]),100);
    [psignal,~]=ksdensity(Feature(signalindx),supportfeat);
    [pnoise,~] =ksdensity(Feature(noiseindx),supportfeat);
    diffpdfs=psignal-pnoise;
    N=length(diffpdfs);
    while diffpdfs(N)>0
        N=N-1;
    end
    if N<length(diffpdfs)
        Th_pdf=supportfeat(N+1);
    else
        Th_pdf=supportfeat(N);
    end
else % if there only SIGNALS
    Th_pdf=min(Feature);
    disp('***************************');
    disp('*                         *');
    disp('*     GREAT VIDEO         *');
    disp('*                         *');
    disp('***************************');
end