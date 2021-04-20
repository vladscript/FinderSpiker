function [SIGNALS,r,H,W,RADroi]=FluorescenceSignals(NC,NumberofVideos,FN,PathName,XY,r)
RADroi={};
%% Load Data *********************************************************
% For each CONDITION and VIDEO
for i=1:NC
    for j=1:str2double(NumberofVideos{i})
        FileName=FN{j,i};
        [mov]=Video_Load(FileName,PathName);    % Load Video
        [FS]=Fluorescence_Load(mov,XY,r);       % Load Fluorescence Signals
        SIGNALS{j,i}=FS;
    end
    disp('***')
end
% Transform to radius cell of ROIs
if iscell(r)
    fprintf('>>Calculating radius of ROIs:')
    RADroi=r;
    RAD=zeros(numel(r),1);
    for xi=1:numel(r)
        PixelROIs=RADroi{xi};
        Xradius=round((max(PixelROIs(:,1))-min(PixelROIs(:,1)))/2);
        Yradius=round((max(PixelROIs(:,2))-min(PixelROIs(:,2)))/2);
        MaxRadius=max([Xradius,Yradius]);
        RAD(xi)=MaxRadius;
    end
    r=RAD;
    fprintf('done.\n')
end
[H,W]=size(mov(1).cdata);   % Height & Width
% clear mov;                  % Clear Video Structure