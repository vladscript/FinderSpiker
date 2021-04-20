function save_activitydata(FolderNamePD,New_Index,Raster_Condition,RASTER_WHOLE_Clean,XY_clean)
% Call Global Variables
global Experiment;
% Direcotry to Save: Up from  this one (pwd)
FileDirSave=getDir2Save();
save([FileDirSave,FolderNamePD,Experiment,'.mat'],'New_Index','Raster_Condition',...
    'RASTER_WHOLE_Clean','XY_clean','-append');
fprintf('>>Automatic Matrix Activity Data from %s: SAVED \n',Experiment)