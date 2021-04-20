%% Script to Load Deafult Directories
% EDIT TO YOUR EXPERIMENTS DIRECTORY:
DefaultPath='C:\Users\Vladimir\Documents\Doctorado\Experimentos\';
if exist(DefaultPath,'dir')==0
    DefaultPath=pwd;
end
%% DO NOT EDIT (preferably):
% Folder with all the files with data and processings: EXPID.mat
FolderNamePD='Processed Data\'; 
% Folder with processings summaries:
FolderNameRT='Resume Tables\';
% Folder with Activity Matrix Features:
FolderNameRaster='Raster Features\';
% Folder with Neural Ensembles Features:
FolderNameEnsembles='Ensemble Features\';
% Folder with Signal Detection Features:
FolderNameDetection='Features Tables\';
% Folder with Network Edges & Nodes Tables:
FolderNameNetwork='NetWorks-CSV\';
% Folder of Datasets:
FolderNameDataset='Experimental Databases\';