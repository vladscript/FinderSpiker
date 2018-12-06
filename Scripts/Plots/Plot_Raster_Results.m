%% PLOT RESULTS
if exist('RASTER_Selected_Clean')
    disp('>>Showing Cleaned Raster')
    Plot_Raster_V(RASTER_WHOLE_Clean,fs);
    RC=R_Condition;
    disp('>>Done.')
else
    disp('>>Showing <preprocessed> Raster')
    Plot_Raster_V(RASTER_WHOLE_Clean,fs);  
    disp('>>Done') % Clean Whole Raster
    RC=Raster_Condition;
end
drawnow;
disp('>>Setting Condition Labels...')
hfig=gcf;
set(hfig,'Name',['ID: ',Experiment(2:end)],'NumberTitle','off')
Label_Condition_Raster(Names_Conditions,RC,fs);   % Labels
disp('>>Ready!')