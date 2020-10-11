%% Setup
Import_FinderSpiker;
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
    % Aski if plots Rate of Activity (RoA)
    okbutton = questdlg('Plot Rate of Activity ?');
    waitfor(okbutton); 
    if strcmp('Yes',okbutton)
        % Plot RoA ********************************************************
        plot_activityrate(R_Condition,Names_Conditions,fs);
        % *****************************************************************
    end
    if numel(R_Condition)>1
        % Aski if plots Rate of Activity (RoA)
        okbutton = questdlg('Plot Subgroups of Neurons by Change of Activity?');
        waitfor(okbutton); 
        if strcmp('Yes',okbutton)
            % Plot Delta RoA Subgroups *************************************
            plot_activitychange(R_Condition,Names_Conditions,fs)
            % **************************************************************
        end
    end
    
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
    fprintf('Execute: \n')
    fprintf('>>Detected_Visual_Inspection\n>>Undetected_Visual_Inspection\n')
    
end
% In case there are Neural Ensemble Analysis
if exist('Features_Condition','var')
    disp('>>Neural Ensembles Ready.')
    okbutton = questdlg('Plot Neural Ensmebles & Hebbian Sequence?');
    waitfor(okbutton); 
    if strcmp('Yes',okbutton)
        % Plot Ensembles **************************************************
        Plot_Neural_Ensembles;
        % *****************************************************************
    end
end