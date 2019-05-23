% Function to accumulate Raster Distances
% Active Cells in each Conditions are neural vectors
% Then, it is calculated the distance among conditions

% Input
%   Load experiment manually from different folder: one by one
%   R_Condition:  Selected rasters of the Experiment
% Output variables @ CSV File
%    Raster_Distance:       Distance Among Conditins
%% Setup
% Initial:
clear; close all; clc;
runs=1;             % Runs Counter
EXPS={};            % List Of Experiments

% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end)),'Processed Data'];
% Load File 
[FileName,PathName,MoreFiles] = uigetfile('*.mat',['Experiment .mat file'],...
    'MultiSelect', 'off',CurrentPathOK);
%% Loop to keep loading files
Raster_Distance=[];
Tbuff=table();
while MoreFiles
    % LOAD DATA
    load([PathName,FileName]);
    [~,NC]=size(R_Condition); % N Conditions
    %%  Loop to Accummulate Data
    if NC>1
        % Matrix of Distances among Conditions
        MatDist=distance_rasters(R_Condition);
        for c=1:NC-1
            for d=c+1:NC
                Raster_Distance(runs,c) = MatDist(c,d);
                S1=Names_Conditions{c};
                S2=Names_Conditions{d};
                S1=S1(isstrprop(S1,'alpha'));
                S2=S2(isstrprop(S2,'alpha'));
                Comparisson{1,c}=[S1,'_',S2];
            end
        end
    else
        disp('>>Not for Single Conditin Experiments')
    end
    % Disp Experiments Selected:
    EXPS{runs,1}=Experiment
    CurrentPathOK=PathName;
    runs=runs+1;
    [FileName,PathName,MoreFiles] = uigetfile('*.mat',['Experiment .mat file'],...
    'MultiSelect', 'off',CurrentPathOK);
end
disp('>>end.')
% Make Table:
Tbuff=table(EXPS,Raster_Distance);
for c=1:NC
    if c==1
        Tbuff.Properties.VariableNames{c}='Experiments';
    else
        Tbuff.Properties.VariableNames{c}=Comparisson{c-1};
    end
end
%% Save Stuff to .mat File at Processed Data Folder

okbutton = questdlg('Save data @ CSV File?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Raster_Distances_',TS,'.csv'];
    % Select Destiny
    % Directory (default)
    CurrentPath=pwd;
    Slshes=find(CurrentPath=='\');
    % [CurrentPath(1:Slshes(end)),'Raster Features']
    CurrentPathOK=[CurrentPath(1:Slshes(end)),'Raster Features'];
    
    PathSave=uigetdir(CurrentPathOK);
    writetable(Tbuff,[PathSave,SaveFile],...
    'Delimiter',',','QuoteStrings',true);

    fprintf('>> Data saved @: %s\n',PathSave)
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')