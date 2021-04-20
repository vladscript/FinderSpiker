% Read coordinates
function [XY,r]=Read_ROI_Coordinates(CurrentPath)
SlahesIndxss=find(CurrentPath=='\');
[XYName,XYPathName] = uigetfile({'*.csv';'*.zip'},['All Coordinates (XY) of ',...
    CurrentPath(SlahesIndxss(end-1)+1:SlahesIndxss(end)-1)],CurrentPath);
if ~isnumeric(XYPathName)
    XYFN = char(XYName);
    switch XYName(end-2:end)
        case 'zip'
            fprintf('>>Reading Coordinates from ImageJ:')
            sROI=ReadImageJROI([XYPathName,XYFN]);
            [XY,r]=getpixelsatROI(sROI);
            % r is a cell array of pixels elliptical ROIs:
        case 'csv'
            fprintf('>>Reading Coordinates from ImPatch: ')
            [x,y,r,~]=textread([XYPathName,XYFN],'%d %d %d %s',...
            'delimiter',',','headerlines',4);
            % r is a vector
            XY = [x,y];
    end
else
    fprintf('No Coordinates Added.\n');
    XY=[];
    r=[];
end
fprintf('done.\n')