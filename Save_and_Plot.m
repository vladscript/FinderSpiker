Import_FinderSpiker;
OddsMatrix=getoddmatrix(isSIGNALS,notSIGNALS);
Experiment=Experiment(Experiment~='\');     % NAMES PATCH

%% Plot Final Data
% GET Condition's Rasters & Whole RASTER 
[~,Raster_Condition,RASTER_CONCAT]=SortNeuronsCondition(RASTER);
% Actual Active ROIs
TotalActiveNeurons=find(sum(RASTER_CONCAT,2)>0);
QoE=round(100*length(TotalActiveNeurons)/length(XY),2);
fprintf('Actual \nActive \nNeurons: %d %%\n',round(QoE));
% SEE RESULTS ################################################
Plot_Raster_Ensembles(RASTER_CONCAT,fs);                           % Clean Whole Raster
% set(gcf,'Name',['ID: ',Experiment],'NumberTitle','off')
CurrentFig=gcf;
CurrentFig.Name = ['ID: ',Experiment,' ready to Analyze'];
CurrentFig.NumberTitle='off';

Label_Condition_Raster(Names_Conditions,Raster_Condition,fs);   % Labels    
%% Update Data: VisualInspector,Oddsmatrix, Raster_Condition
checkname=1; oksave=false;
while checkname==1
    % Get Directory
    DP=pwd;
    Slashes=find(DP=='\');
    DefaultPath=[DP(1:Slashes(end)),'Processed Data'];
    if exist(DefaultPath,'dir')==0
        DefaultPath=pwd; % Current Diretory of MATLAB
    end
    uiwait(msgbox('Open .mat file to overwrite & update changes'))
    [FileName,PathName] = uigetfile('*.mat',[' SAVE CHANGES & pick the Analysis File ',Experiment],...
        'MultiSelect', 'off',DefaultPath);
    % dotindex=find(FileName=='.');
    if strcmp(FileName(1:end-4),Experiment)
        checkname=0; oksave=true;
        % SAVE DATA
        save([PathName,FileName],'OddsMatrix','Raster_Condition',...
            'RasterAlgorithm','RASTER_CONCAT',...
            '-append');
        disp([Experiment,'   -> UPDATED (OddsMAtrix Data)'])
    elseif FileName==0
        checkname=0;
        disp('....CANCELLED')
    else
        disp('Not the same Experiment!')
        disp('Try again!')
    end
end    
%% Save Features Tables
if oksave
    save_preocessing_intel;
end