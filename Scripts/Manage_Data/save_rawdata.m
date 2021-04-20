function save_rawdata(FolderNamePD,NumberofVideos,XY,r,PathName,RADroi,dyename,FN)
% Call Global Variables
global SIGNALS;
global Names_Conditions;
global fs;
global Experiment;
% Directory to Save: Up from  this one (pwd)
FileDirSave=getDir2Save();
% Creata Folder:
if ~isdir([FileDirSave,FolderNamePD])
    fprintf('Folder [%s] created @ %s\n',FolderNamePD,FileDirSave)
    mkdir([FileDirSave,FolderNamePD]);
end
% Saves data (and the day)
save([FileDirSave,FolderNamePD,'\',Experiment,'.mat'],'Experiment','SIGNALS',...
    'Names_Conditions','NumberofVideos','XY','fs','r','PathName','RADroi',...
    'dyename','FN');
fprintf('>>Raw Fluorescence Data from %s: SAVED \n',Experiment)