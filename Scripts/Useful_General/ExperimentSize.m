function [fs,dyename,NV,NC]=ExperimentSize(NumberofVideos,FN)
%% Read Sampling Frequency
fs=NaN; % To read fs to get into the following loop:
NV=1;
while isnan(fs) % Interface error user reading fs
    fs = inputdlg('Sampling Frequency [Hz] : ',...
                 'Frames per Second or', [1 75]);
    fs = str2double(fs{:});
end
% Read Fluorophore DYe
dyename = inputdlg('Fluorophore : ',...
             'DYE', [1 75]);

for v=1:length(NumberofVideos)
    NVal(v)=round(str2double(NumberofVideos(v)));
end
NV=max(NVal);   % Max N of Videos
NC=size(FN,2);                                % N Conditions