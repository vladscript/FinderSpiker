%% Load Datasetes ******************************************************

Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end))]; % FinderSpiker Root Folder

%% Raster Activity Minimal Constraint:
[FileName,PathName] = uigetfile('*.csv','Raster Activity Dataset',...
    'MultiSelect', 'off',CurrentPathOK);
CurrentPathOK=PathName;
Xraster=readtable([PathName,FileName]);
Yraster=categorical(table2array( Xraster(:,1)) );   % Labels
EXPIDraster=table2array( Xraster(:,2));             % EXPIDs Cell of Strings
FeaturesRaster=Xraster.Properties.VariableNames;    % Feature Names

%% Ensembles General Minimal Constraint:
[FileName,PathName] = uigetfile('*.csv','General Ensembles Dataset',...
    'MultiSelect', 'off',CurrentPathOK);
if FileName==PathName
    Xensemble=table;        % empty table
    FeaturesEnsemble={};    % empty cell
    EXPIDensemble={};       % empty cell
    Yensemble=[];           % empty vector
    disp('>>No Ensemble Features');
else
    disp('>>Reading Ensemble Features');
    Xensemble=readtable([PathName,FileName]);
    Yensemble=categorical(table2array( Xensemble(:,1)) );   % Categorical ARRAYS
    EXPIDensemble=table2array( Xensemble(:,2));             % EXPIDs Cell of Strings
    FeaturesEnsemble=Xensemble.Properties.VariableNames;    % Feature Names
    ignftr=[15,16,17]; % Links/Neuron A->euron B Most Connected
    acceptcols=setdiff(1:numel(FeaturesEnsemble),ignftr);
    Xensemble=Xensemble(:,acceptcols);
    FeaturesEnsemble=FeaturesEnsemble(acceptcols);
    disp('>>Done.');
end



%% Load NETWORK csv
[FileName,PathName,MoreFiles] = uigetfile('*.csv','Network Dataset file',...
    'MultiSelect', 'off',CurrentPathOK);
Xnet=readtable([PathName,FileName]);
Ynet=categorical(Xnet{:,1} );     % Labels
EXPIDnet=Xnet{:,2};               % Cell of Strings
FeaturesNet=Xnet.Properties.VariableNames;      % Feature Names

%% MERGE DATASETS: RASTER|ENSEMBLES|NETWORK
% get same sorting for labels and Dateset rows 
[YlabelsA,indxRasterA,indxEnsembles]=getindxofAinB([EXPIDraster,Yraster],[EXPIDensemble,Yensemble]);
[YlabelsB,indxRasterB,indxNet]=getindxofAinB([EXPIDraster,Yraster],[EXPIDnet,Ynet]);
% Check New Indexing
C=intersect(YlabelsA,YlabelsB,'rows');
if size(YlabelsA,1)==size(C,1) && size(YlabelsA,2)==size(C,2)
    disp('>>Indexing: OK')
    Ylabels=YlabelsA;

    % MERGE TABLE
    X=[table(Ylabels(:,2),Ylabels(:,1)),Xraster(indxRasterA,3:end),...
        Xensemble(indxEnsembles,3:end),...
        Xnet(indxNet,3:end)];
    X.Properties.VariableNames{1}='Condition';
    X.Properties.VariableNames{2}='EXP_ID';
    % FILTER REJECTED EXPERIMENTS
    NData=size(X,1);

    %% Make and Save Table
    okbutton = questdlg('Make CSV Table?');
    waitfor(okbutton); 
    if strcmp('Yes',okbutton)
        % Set Save Name
        CurrentPathOK=pwd;
        timesave=clock;
        TS=num2str(timesave(1:5));
        TS=TS(TS~=' ');
        SaveFile=['\ALL_Features_DATASET_',TS,'.csv'];
        % Select Destiny
        PathSave=uigetdir(CurrentPathOK,'Set Directory to save ALL Features Table');
        disp('>>Making CSV table...')
        writetable(X,[PathSave,SaveFile],...
                        'Delimiter',',','QuoteStrings',true);
        fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
    else
        fprintf('>>Unsaved data.\n')
    end
else
    disp('>>ERROR: Check Tables')
end
% WITHOUT THE REJECTED ONES