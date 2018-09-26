% Functio to calculate time durations of (selected) raster conditions
% Input
%   Onsets (frames)
%   R_Condition (Cells x Frames) Cell of Nc Conditions
%   fs: Sampling Frequency (Hz)
% Output
%   Table of Durations
function RasterDurations=get_raster_durations(Onsets,R_Condition,fs)
% Setup
Nc=numel(R_Condition);
for c=1:Nc
    [~,Framces]=size(R_Condition{c});
    InitialMinute=(Onsets{c}-1)/fs/60;
    RasterDurations(c,:)=[InitialMinute,InitialMinute+Framces/fs/60];
end