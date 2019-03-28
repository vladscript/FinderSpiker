% Function that gets the pdf aprox from the histogram of 
% Inter (Ca++) Transient Interval histogram obtained from the raster
% Input
%   R: raster data: Cells x Frames
%   fs: sampling frequency
% Ouput
%   ISIbin: bins of ISI [seconds]
%   ISIp:   probaility of ISI
%   Lengthbin: bins of Transient Length [seconds]
%   Lengthp:   probaility of Transient Length
%   StatsFeatures: Vector of Descriptive Statistics:
%   MEAN, MODE, VARIANCE,SKEWNESS,KURTOSIS
%   for ITI & TrasnientLEngth 
function [ISIbin,ISIp,Lengthp,Lengthbin,StatsFeatures]=get_iti_pdf(R,varargin)
%% Setup
[Cells,~]=size(R);
if isempty(varargin)
    disp('No Sampling Frequency->Discrete Time Domain')
    fs=1;
else
    fs=varargin{1};
end
%% Main Loop
ISIs=[];
TranLengths=[];
for c=1:Cells
    r=R(c,:);
    [ISIcell,TranLengthscell]=interval_duration_events(r);
    ISIs=[ISIs,ISIcell];
    TranLengths=[TranLengths,TranLengthscell];
    % Transien Rate
    RoT(c)=numel(TranLengthscell)*fs/numel(r)*60; % Ca Transients per MINUTE
    disp(c);
end
%% GET PDFs
if numel(ISIs(ISIs>0))>1
    [ISIp,ISIbin]=ksdensity(ISIs/fs,'support','positive','function','survivor');
else
    ISIbin=[];
    ISIp=[];
end
if numel(TranLengths(TranLengths>0))>1
    [Lengthp,Lengthbin]=ksdensity(TranLengths/fs,'support','positive','function','survivor');
else
    Lengthbin=[];
    Lengthp=[];
end
    
%% STATISTICS
%   MEAN, MODE, MEDIAN, VARIANCE,SKEWNESS,KURTOSIS x ITI and Length pdfs
StatsFeatures=[mean(ISIs/fs),mode(ISIs/fs),median(ISIs/fs),var(ISIs/fs),skewness(ISIs/fs),kurtosis(ISIs/fs),...
mean(TranLengths/fs),mode(TranLengths/fs),median(TranLengths/fs),var(TranLengths/fs),skewness(TranLengths/fs),kurtosis(TranLengths/fs),...
mean(RoT),mode(RoT),median(RoT),var(RoT),skewness(RoT),kurtosis(RoT)];
