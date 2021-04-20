function [SIGNALS,r,H,W,RADroi,fs]=FluorescenceSignalsCSV(NumberofSplits,FN,PathName,r)
RADroi={};
%% Load Data *********************************************************
% For each CONDITION and VIDEO
[X,fs]=readdeltaFzero([PathName,FN{1}]);
[C,F]=size(X);
if C>F
    X=X';
    [C,F]=size(X);
    disp('Transposed to Cells x Frames');
end
FramesSeg=round(linspace(1,F+1,NumberofSplits+1));

for i=1:1
    for j=1:NumberofSplits
        fprintf('>>Splitting Data %i of %i\n',j,NumberofSplits);
        SIGNALS{j,i}=X(:,FramesSeg(j):FramesSeg(j+1)-1);
    end
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

name='Frame Size Pixels';
numlines=[1 75];
Names_Conditions=inputdlg({'Height';'Weight'},name,numlines,{'100';'200'});
H=str2double( Names_Conditions{1} );
W=str2double( Names_Conditions{1} );
% [H,W]=size(mov(1).cdata);   % Height & Width
% clear mov;                  % Clear Video Structure