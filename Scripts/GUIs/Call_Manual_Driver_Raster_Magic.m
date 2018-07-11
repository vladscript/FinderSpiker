%% OLD STUFF
button = questdlg('Results Inspection?');
if strcmp('Yes',button)    
    [NV,NC]=size(isSIGNALS);
    preisSIGNALS=isSIGNALS;
    [isSIGNALSOK,SIGNALSclean,DRIVEROK,RASTEROK,...
    LAMBDASSpro,SNRlambda,OddsMatrix]=Manual_Driver_Raster_Magic(isSIGNALS,SIGNALS,...
    DETSIGNALS,preDRIVE,preLAMBDAS,RASTER,Responses,Names_Conditions,fs);
    
    % Re-Sort WHOLE RASTER 
    [New_Index,Raster_Condition,RASTER_WHOLE]=SortNeuronsCondition(RASTEROK);
    RASTER_WHOLE_Clean=RASTER_WHOLE(New_Index,:);
    XY_clean=XY(New_Index,:);
    % Clean Raster and Coordinates
    TotalActiveNeurons=find(sum(RASTER_WHOLE_Clean,2)>0);                % INDEX of Active NEURONS
    QoE=round(100*length(TotalActiveNeurons)/length(XY),2);
    fprintf('Actual Active Neurons: %d %%\n',round(QoE));
    % Whole Raster
    RASTER_WHOLE_Clean=RASTER_WHOLE_Clean(TotalActiveNeurons,:);
    XY_clean=XY_clean(TotalActiveNeurons,:);                        % Clean Coordinates
    New_Index_Active=New_Index(TotalActiveNeurons);
    % SEE RESULTS ################################################
    Plot_Raster_V(RASTER_WHOLE_Clean,fs);                           % Clean Whole Raster
    set(gcf,'Name',['ID: ',Experiment(2:end)],'NumberTitle','off')
    Label_Condition_Raster(Names_Conditions,Raster_Condition,fs);   % Labels    
    % Update Results                [ok] ----------------------------------
    % Ask for Directory to save & MAT file to update
    checkname=1;
    while checkname==1
        DefaultPath='C:\Users\Vladimir\Documents\Doctorado\Software\GetTransitum\Calcium Imaging Signal Processing\FinderSpiker\Processed Data';
            if exist(DefaultPath,'dir')==0
                DefaultPath=pwd; % Current Diretory of MATLAB
            end

        [FileName,PathNamePro] = uigetfile('*.mat',[' Pick the Analysis File ',Experiment],...
            'MultiSelect', 'off',DefaultPath);
        dotindex=find(FileName=='.');
        if strcmp(FileName(1:dotindex-1),Experiment(2:end))
            checkname=0;
            % SAVE DATA
            save([PathNamePro,FileName],'isSIGNALSOK','SIGNALSclean','DRIVEROK','SNRlambda','LAMBDASSpro',...
                'New_Index','Raster_Condition','RASTER_WHOLE_Clean','XY_clean','RASTEROK','OddsMatrix','-append');
            disp([Experiment,'   -> UPDATED: Visual Review'])

        else
            disp('Not the same Experiment!')
            disp('Try again!')
        end
    end   
end