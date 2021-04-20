% Retrieve Single Directory
% Input: number; output: sufix for:
% 1: Processed
% 2: Resume
% 3: Raster Features
% 4: Ensembles Features
% 5: Signal Features Tables
% 6: NetworksCSVs
% 7: Experimental Dataset
function DirectorySufix=DirectorySelector(numDir)
DirectorySufix=pwd;
Load_Default_Directories;
switch numDir
    case 1
        DirectorySufix=FolderNamePD;
    case 2
        DirectorySufix=FolderNameRT;
    case 3
        DirectorySufix=FolderNameRaster;
    case 4
        DirectorySufix=FolderNameEnsembles;
    case 5
        DirectorySufix=FolderNameDetection;
    case 6
        DirectorySufix=FolderNameNetwork;
    case 7
        DirectorySufix=FolderNameDataset;
    otherwise
        disp('No directory selected')
end