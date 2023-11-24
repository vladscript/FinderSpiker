%% Make Data set of Statistics of Network Features
%% 0. Load MAT File
fprintf('\n>Read node parameters from .mat file created with Get_Gephi_Data\n')
fprintf('to compute numerical statistics for condition:\n mean\n median\n mode\n variance\n kurtosis\n skewness\n')
fprintf('It considers Nodes with Degree>0 in C_i or peviously Degree>0 in C_i-1\n');
Load_Default_Directories;
% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end)),FolderNameNetwork];
[FileName,PathName] = uigetfile('.mat','Load MATLAB file of GEPHI Node Features',...
    CurrentPathOK);
fullFileName = fullfile(PathName, FileName);
load(fullFileName);
% setup
% [NE,NC]=size(GEPHIDATA); % Number of Experiments & Conditions
MoreFeats=true; % enter the loop

%% 1. Chooose Feature from table: clustering, grade, etc, etc
busca=true;
auxn=1;
while busca
    if ~isempty( GEPHIDATA{auxn} )
        NetFeaturesNames=GEPHIDATA{1}.Properties.VariableNames(3:end)';
        busca=false;
    else
        auxn=auxn+1;
    end
end
while MoreFeats    
    [index_var,index_CHECK] = listdlg('PromptString','Select a Network Feature:',...
                        'SelectionMode','single',...
                        'ListString',NetFeaturesNames);
    ActualFeature=NetFeaturesNames(index_var);
    %% 1 Check Empty Tables and Empty-Feature Columns

    [ExsitFeature,EXPLIST]=emptycolnettables(GEPHIDATA,Names_Conditions,ActualFeature);
    %% 2. Read data from ALL Tables and ALL Conditions

    [TotalTable,DATAnet]=TableNetwork(ExsitFeature,GEPHIDATA,EXPLIST,Names_Conditions,ActualFeature);
    

    %% 3. Make Table for ML o Statistical Anlysis
    % Set Name to Variables
    TotalTable.Properties.VariableNames={'Condition','EXPID',...
        ['mean_',ActualFeature{1}],['mode_',ActualFeature{1}],...
        ['median_',ActualFeature{1}],['var_',ActualFeature{1}],...
        ['skew_',ActualFeature{1}],['kurt_',ActualFeature{1}]};

    
    okbutton = questdlg('Make CSV Table?');
    waitfor(okbutton); 
    if strcmp('Yes',okbutton)
        % Netwrok Directory
        CurrentPath=pwd;
        Slshes=find(CurrentPath=='\');
        % CurrentPathOK=[Dirpwd(1:slashesindx(end)),'NetWorks-CSV'];
        CurrentPathOK=[Dirpwd(1:slashesindx(end)),FolderNameNetwork];
        % Set Save Name
        timesave=clock;
        TS=num2str(timesave(1:5));
        TS=TS(TS~=' ');
        SaveFile=['Network_',ActualFeature{1},'_Dataset_',TS,'.csv'];
        % Select Destiny
        PathSave=uigetdir(CurrentPathOK,'Choose destination folder:');
        disp('>>Making CSV table...')
        writetable(TotalTable,[PathSave,SaveFile],...
                        'Delimiter',',','QuoteStrings',true);
        % fprintf('>> Dataset saved @: %s\n',[PathSave,SaveFile])
        fprintf('<a href="matlab:dos(''explorer.exe /e, %s, &'')">See network parameters here</a>\n',PathSave);
    else
        fprintf('>>Unsaved dataset.\n')
    end

    answer = questdlg('Inspect More Features?', ...
    'NETWROK FEATURES', ...
    'YES','NO','NO');

    switch answer
        case 'YES'
            disp([answer ' more Network features coming right up.'])
            MoreFeats = true;
        case 'NO'
            disp([answer ' more Network features.'])
            MoreFeats = false;
    end

    % Ask if Save CSV, clear and Goodbye
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n Run next: \n>Concatenate_NetFeats\n to concatenate tables of experiments\n');
%% END OF THE WORLD