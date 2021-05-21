% Read coordinates
function [XY,r]=Read_ROI_Coordinates(CurrentPath)
DeafultROIRadius=5;
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
            
            fileID = fopen([XYPathName,XYFN],'r');
            fprintf('\n> Reading coordinates from: ')
            IsNumberRow2 = textscan(fileID, '%f%f',1, 'Delimiter',...
                ',', 'HeaderLines', 1, 'ReturnOnError', true);
            if isempty(IsNumberRow2{1})
                fprintf('ImPatch\n')
                [x,y,r,~]=textread([XYPathName,XYFN],'%d %d %d %s',...
                'delimiter',',','headerlines',4);
                % r is a vector
                XY = [x,y];
            else
                fprintf('Row-2 CSV UTF 8 with Byte Order Mark\n')
                fclose(fileID);
                % Read from coordinates @row 2 CSV
                startRow=2;
                XY = readXYfromCSV([XYPathName,XYFN],startRow);
                fprintf('Setting default ROI radius to %i pixels \n',DeafultROIRadius)
                r=DeafultROIRadius*ones(size(XY,1),1);
            end
            
    end
else
    fprintf('No Coordinates Added.\n');
    XY=[];
    r=[];
end
fprintf('done.\n')