% Find Peaks
% Find peaks from signal by a given threshold.
%
% [PeaksIdx] = FindPeaks_JP(X,th,join)
%
% Inputs
% X = data as vector Fx1 (F = #frames)
% th = threshold
% join = set mode to get peaks (0 = each vector above threshold is a peak;
%        1 = joining of adjacent vectors above the threshold is a peak)
% 
% Outputs
% PeaksIdx = Fx1 vector containing the peak indices
%
% ..:: by Jesús E. Pérez-Ortega ::.. Feb-2012
% JP debug 1-oct-2012
% JP add 3rd input to join or not

function [PeaksIdx] = FindPeaks_JP(X,th,join)

F = numel(X);
PeaksIdx=zeros(F,1);

% Get index
idx = find(X>=th);                         % index above threshold
peaks=numel(idx);                          % peaks

if ~peaks
    disp('No data found')
    return
end
for i=1:peaks
    PeaksIdx(idx(i))=i;
end

% Join peaks
if join
    is = find(idx~=[0; idx(1:numel(idx)-1)+1]);    % index of same peak

    % Delete first if start above threshold
    if min(idx)==1
        if numel(is)>1
            idx=idx(is(2):numel(idx));
            PeaksIdx(1:is(2)-1)=0;
            is=is(2:numel(is))-is(2)+1;
        else
            disp('No data found')
            return
        end
    end

    % Delete last if ends above threshold
    if(max(idx)==F)
        if numel(is)>1
            idx=idx(1:max(is)-1);
            PeaksIdx(is(numel(is)-1):F)=0;
            is=is(1:numel(is)-1);
        else
            disp('No data found')
            return
        end
    end

    % number of total peaks
    peaks = numel(is);                                       
    if peaks
        for j = 1:peaks-1
            PeaksIdx(idx(is(j)):idx(is(j+1)-1),1)=j;    % set #peak
        end
        PeaksIdx(idx(is(peaks)):max(idx),1)=peaks;
    end
end