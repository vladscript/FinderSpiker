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
    [ISIs,TranLengths]=interval_duration_events(r);
    % disp(c);
%     if ~isempty(r(r>0))
%         dr=diff(r); 
%         % +dr -> start Interval
%         % -dr -> end Interval
%         StartsIntervasl=find(dr>0);
%         EndsIntervasl=find(dr<0);
%         if and(~isempty(StartsIntervasl),~isempty(EndsIntervasl))
%             % If it started at Initial Transient
%             if StartsIntervasl(1)>EndsIntervasl(1)
%                 EndsIntervasl=EndsIntervasl(2:end);
%             end
%             % If it ended at with out End of Transient
%             if StartsIntervasl(end)>EndsIntervasl(end)
%                 StartsIntervasl=StartsIntervasl(1:end-1);
%             end
%             if numel(StartsIntervasl)==numel(EndsIntervasl)
%                 if numel(StartsIntervasl)>1
%                     ISIs=[ISIs,StartsIntervasl(2:end)-EndsIntervasl(1:end-1)];
%                 end
%                 TranLengths=[TranLengths, EndsIntervasl-StartsIntervasl];
%             else
%                 disp('WTF?');
%             end
%         end
%     end
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
%   MEAN, MODE, VARIANCE,SKEWNESS,KURTOSIS x ITI and Length pdfs
StatsFeatures=[mean(ISIs/fs),mode(ISIs/fs),var(ISIs/fs),skewness(ISIs/fs),kurtosis(ISIs/fs),...
mean(TranLengths/fs),mode(TranLengths/fs),var(TranLengths/fs),skewness(TranLengths/fs),kurtosis(TranLengths/fs)];
% %% CHECK PLOT
% figure
% subplot(1,2,1); plot(ISIbin,ISIp)
% axis tight; grid on;
% title('InterTranInterval pdf')
% xlabel('[s]')
% subplot(1,2,2); plot(Lengthbin,Lengthp)
% axis tight; grid on;
% title('Ca++ Tran. Length pdf')
% xlabel('[s]')