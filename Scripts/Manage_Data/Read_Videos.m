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
%   r: Radious of Coordinates

function [Names_Conditions,NumberofVideos,FN,PathName,XY,r]=Read_Videos(DefaultPath)
%% Ask how many conditions and how many videos *********************
% 1st Input DIalogue
NC = inputdlg('Number of Conditions: ',...
             'VIDEOS', [1 25]);
NC = str2num(NC{:}); 
Conditions_String='Condition_';
n_conditions=[1:NC]';
Conditions=[repmat(Conditions_String,NC,1),num2str(n_conditions)];
Cond_Names=cell(NC,1);
% Names=cell(NC,1);
for i=1:NC
    Cond_Names{i}=Conditions(i,:);
    Names_default{i}=['...'];
    NVids_default{i}=['2'];
end
% 2nd Input Dialogue
name='Names';
numlines=[1 25];
Names_Conditions=inputdlg(Cond_Names,name,numlines,Names_default);
% 3th Input Dialogue
name='Videos';
numlines=1;
NumberofVideos=inputdlg(Names_Conditions,name,numlines,NVids_default);

%% Read and get all Coordinates*********************************************
CurrentPath=DefaultPath;

for i=1:NC
    % Read Videos
    [FileName,PathName] = uigetfile('*.avi',[NumberofVideos{i},' VIDEOS ',Names_Conditions{i}],...
    'MultiSelect', 'on',CurrentPath);
    FileName=char(FileName);
    for j=1:str2double(NumberofVideos{i})
        FN{j,i} = FileName(j,:); % Different Names of Videos
    end
    CurrentPath=PathName;
end
% Read coordinates
[XYName,XYPathName] = uigetfile('*.csv',['All Coordinates (XY) ',Names_Conditions{i}],CurrentPath);
XYFN = char(XYName);
[x,y,r,~]=textread([XYPathName,XYFN],'%d %d %d %s','delimiter',',','headerlines',4);
XY = [x,y];