%% Setup
Update_Directory;
Experiment=Experiment(Experiment~='\');     % NAMES PATCH

%% Script that plots the raster of the Selected OR the Whole Clean raster
if exist('RASTER_Selected_Clean','var')
    disp('>>Plotting Selected Raster to Analyze: ')
    Plot_Raster_Ensembles(RASTER_Selected_Clean,fs);
    CurrentFig=gcf;
    CurrentFig.Name=['ID: ',Experiment,' Selected to Analyze'];
    CurrentFig.NumberTitle='off';
    
    Label_Condition_Raster(Names_Conditions,R_Condition,fs);   % Labels
    disp('>>Raster Ready to Analyze.')
elseif exist('RASTER_CONCAT','var')
    disp('>>Plotting Clean Raster: ')
    Plot_Raster_Ensembles(RASTER_CONCAT,fs);
    CurrentFig=gcf;
    CurrentFig.Name=['ID: ',Experiment,' Clean Raster'];
    CurrentFig.NumberTitle='off';
    Label_Condition_Raster(Names_Conditions,Raster_Condition,fs);   % Labels
    disp('>>Next Step: Select Frames of Raster.')
else 
    disp('>>Plotting Raster from Automatic Processing ')
    Plot_Raster_Ensembles(RASTER_WHOLE_Clean,fs);
    CurrentFig=gcf;
    CurrentFig.Name=['ID: ',Experiment,' from Automatic Processing'];
    CurrentFig.NumberTitle='off';
    Label_Condition_Raster(Names_Conditions,RASTER,fs);   % Labels
    disp('>>Next Step: Inspect Signals: Detected & Undetected.')
end
