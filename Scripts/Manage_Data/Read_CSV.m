% Funtion that reads Names and Directory of Videos
% and loads coordinates
% Input (same directory) from dialogue: 
%   ALL (XY).csv
%   Condition_A_00.avi
%   Condition_A_01.avi
%   ...
%   Condition_B_00.avi
%   ... 
%   Condition_Z_XX.avi
% Output
%   Names_Conditions
%   NumberofVideos (per condition)
%   FN{}: Structure with the names of the videos
%   PathName: Directory of the videos
%   XY: Coordinates of Active Cells
% If ImPatch ROI setlist
%       r: Radious of Coordinates  
% If ImageJ ROI setlist
%       r: Cell Array of Coordinates in ROI
function [Names_Conditions,NumberofVideos,FN,PathName,XY,r]=Read_CSV(DefaultPath)
%% Ask how many conditions and how many videos *********************
NC = 1; 
Conditions_String='Condition_';
n_conditions=[1:NC]';
Conditions=[repmat(Conditions_String,NC,1),num2str(n_conditions)];
Cond_Names=cell(NC,1);
% Names=cell(NC,1);
for i=1:NC
    Cond_Names{i}=Conditions(i,:);
    Names_default{i}=['...'];
    NVids_default{i}=['1'];
end
% 2nd Input Dialogue
name='Condition Name';
numlines=[1 85];
Names_Conditions=inputdlg(Cond_Names,name,numlines,Names_default);
% 3th Input Dialogue
name=' Split CSV data in N Segments:';
% numlines=1;
NumberofVideos=inputdlg(Names_Conditions,name,numlines,NVids_default);

%% Read and get all Coordinates*********************************************
% CurrentPath=DefaultPath;

% Read Videos
[FN{1},PathName] = uigetfile('*.csv','Select CSV file with Fluorescence Traces',...
    'MultiSelect','off',DefaultPath);
[XY,r]=Read_ROI_Coordinates(PathName);

% % Read coordinates
% [XYName,XYPathName] = uigetfile({'*.csv';'*.zip'},['All Coordinates (XY) ',Names_Conditions{i}],CurrentPath);
% if ~isnumeric(XYPathName)
%     XYFN = char(XYName);
%     switch XYName(end-2:end)
%         case 'zip'
%             fprintf('>>Reading Coordinates from ImageJ:')
%             sROI=ReadImageJROI([XYPathName,XYFN]);
%             [XY,r]=getpixelsatROI(sROI);
%             % r is a cell array of pixels elliptical ROIs:
%         case 'csv'
%             fprintf('>>Reading Coordinates from ImPatch: ')
%             [x,y,r,~]=textread([XYPathName,XYFN],'%d %d %d %s','delimiter',',','headerlines',4);
%             % r is a vector
%             XY = [x,y];
%     end
% else
%     fprintf('No Coordinates Added.\n');
%     XY=[];
%     r=[];
% end
% fprintf('done.\n')