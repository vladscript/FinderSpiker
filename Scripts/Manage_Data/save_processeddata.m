function save_processeddata(FolderNamePD,ESTSIGNALS,SNRwavelet,TAUSall,isSIGNALS,notSIGNALS)
% Call Global Variables:
global DETSIGNALS;
global preDRIVE;
global preLAMBDAS;
global RASTER
global Responses
global SIGNALSclean;
global SNRlambda;
global Experiment;
global RasterAlgorithm;
% Direcotry to Save: Up from  this one (pwd)
FileDirSave=getDir2Save();
save([FileDirSave,FolderNamePD,Experiment,'.mat'],'DETSIGNALS','ESTSIGNALS',...
    'SNRwavelet','SIGNALSclean','SNRlambda','RasterAlgorithm',...
    'preDRIVE','preLAMBDAS','TAUSall','RASTER','isSIGNALS','notSIGNALS',...
    'Responses','-append');
fprintf('>>Processed Fluorescence Data from %s: SAVED \n',Experiment)