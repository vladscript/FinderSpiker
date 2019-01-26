% Script to sort  and merge for Rasters for Several Analized Conditions
% from get_bayes_ensembles and from NeuralNetworks by JP
% Input
%  NC:  How many Conditions: Input Dialogue
%  Raster_Analays Strucutre from NN: [ Frames x Cells ]
%  Strucutre from NN:   input from workspace:
%                       Raster
%                       Significative Frames
%                       Label of Ensembles
%  RASTER_Selected_Clean: to get Total SELECTED Frames
%  XY_Selected: Selected Active Coordinates
%  fs
%  Experiment ID
% Output
%  New_Order_Clustering:    Ensemble Sorting
%  XY_cluster:              Ensembled-Sorted Indexes
%  Neural Ensemble Features:
%               Nensembles    
%  Functional Network Features:                
%               Gephi CSV files

%% Setup
Update_Directory;
%% Select Number of Analyzed Rasters ######
NC = inputdlg('Enter Number of Analyzed Conditions:',...
         'Conditions', [1 50]);
NC= str2double(NC{:}); 

%% Read Result from NeuralNetworks Clustering Results Structure
WorkspaceVariables=who;     % Variable Names at Workspace
Condition_Names={};         % Condition Names
% Filter Only with term *[_Analysis]*
Filterby='_Analysis'; % Sufix for Analysis Variables
[index_var,WorkspaceVariables]=get_var_index(WorkspaceVariables,Filterby,NC);
% DATA ANALYSIS AND NAMES CONDITIONS
for c=1:NC
    % Get Structure Data Clustering Results
    VarNames{c}=WorkspaceVariables{index_var(c)}; % USED TO SAVE ANALYSIS VAR
    NN_data_Read{c}=evalin('base',WorkspaceVariables{index_var(c)}); % Data from WORKSPACE
    % Get Conditon Names from Variables
    AuxName=WorkspaceVariables{index_var(c)};              % Get Conditon: R_'Xcondition'
    indexes_name=find(AuxName=='_');
    Condition_Names{c}=AuxName(indexes_name(1)+1:indexes_name(end)-1);
end
%% Get DATA from Analysis Variables and CONCATENATE
% RASTERS ARE: [FRAMES x CELLS] -----------
[Rasters,LENGHTRASTER,NGroups,SIG_FRAMES,LABELS,THR,...
 ExperimentRaster,signif_frames,labels_frames,CummFrames,CummGroups,...
 TotalNG]=get_ensemble_intel(NN_data_Read);
clear NN_data_Read;

%% Check if it was Analyzed the Complete Experiment
TotalFrames=size(RASTER_Selected_Clean,2); % INPUT
IndexesActive = find(sum(ExperimentRaster));    % Active Neurons
if CummFrames==TotalFrames
    disp('Complete Experiment')
    ExperimentRasterClean = ExperimentRaster;
    XY_selectedClean = XY_selected;
else
    % Select Coordinates from Concatenated Raster (if is not complete)
    % Update & Clean Data:
    disp('Partial Experiment')
    ExperimentRasterClean = ExperimentRaster(:,IndexesActive);
    XY_selectedClean = XY_selected(IndexesActive,:);
end

Indexes=1:length(XY_selectedClean); % Indexes for Active Neurons
Index_Ensemble=Indexes;
CAG=sum(ExperimentRasterClean,2);   % CAG :Co-Activity-Graphy
%% Sorting Clustering
[New_Order_Clustering,Neurons_State]=OrderClusters(labels_frames,signif_frames,ExperimentRasterClean,TotalNG);
% Coordinates
XY_cluster=XY_selectedClean(New_Order_Clustering,:); % Re-SORTED COORDINATES OF ENSEMBLES
% Indexes=sort(Indexes(New_Order_Clustering));    % useless
%% COLORMAP ENSEMBLES
ColorState=colormapensembles(TotalNG,NC,NGroups);
%% Plot Ensembles of Whole Raster ---------------------------------------------------------------------
Experiment=Experiment(Experiment~='\');     % NAMES PATCH
%   Original ****************************************
OriginalExperiment=ExperimentRasterClean;
Plot_Raster_Ensembles(OriginalExperiment',fs,5,Indexes);    % Disorted Raster
disp('Coloring Ensembles...')
Plot_State_Colors(labels_frames,signif_frames,ColorState,OriginalExperiment,fs,CAG,Indexes);
disp('Coloring Ensembles Done.')
% plot_CAG_threshold(THR,R_Condition,fs)
plot_CAG_threshold(THR,LENGHTRASTER,fs)
if CummFrames==TotalFrames
    Label_Condition_Raster(Names_Conditions,R_Condition,fs);   % Labels
end
Figg=gcf; Figg.Name=['Neural Ensembles of ',Experiment];
%   Sorted *******************************************
Plot_Raster_Ensembles(OriginalExperiment',fs,1,Indexes(New_Order_Clustering));   % Sorted Raster
% Plot_State_Colors;
disp('Coloring Ensembles...')
Plot_State_Colors(labels_frames,signif_frames,ColorState,OriginalExperiment,fs,CAG,Indexes(New_Order_Clustering));
disp('Coloring Ensembles Done.')
plot_CAG_threshold(THR,LENGHTRASTER,fs)
if CummFrames==TotalFrames
    Label_Condition_Raster(Names_Conditions,R_Condition,fs);   % Labels
end
Figg=gcf; Figg.Name=['Neural Ensembles (resorted) of ',Experiment];
% Ensemble Transitions HEBBIAN SEQUENCE **************
Ensembles_Transitions(fs,labels_frames,signif_frames,ColorState,1);

%% IF ALL EXPERIMENT in ONE Analysis
NCplot=NC;
UniRMutiE=false; % Checker if all conditons are in one single analyzed raster
if CummFrames==TotalFrames && NC==1 && numel(Names_Conditions)>1
    disp('---------------------------------------------------------------')
    disp('>>>>>>>>>>>>>   Extracting Conditons from    <<<<<<<<<<<<<<<<<<')
    disp('>>>>>>>>>>>>>>>>   Whole Experiment   <<<<<<<<<<<<<<<<<<<<<<<<<')
    disp('---------------------------------------------------------------')
    NCplot=length(Names_Conditions);
    THR=repmat(THR,NCplot,1);
    UniRMutiE=true;
    NCupdate=length(Names_Conditions);
    cummFrames=1;
    intA=0; % Initialize Interval A
    intB=0; % Initialize Interval B
    % DISCLOSE IF WHOLE RASTER was ANALYZED
    for c=1:NCupdate
        Rasters{c}=R_Condition{c}';        % frames x cells
        [~,framesC]=size(R_Condition{c});  % frames of condition
        intB=intB+framesC;
        [SIG_FRAMES_UPD{c},indexSIG,~]=intersect(SIG_FRAMES{1},intA+1:intB);
        SIG_FRAMES_UPD{c}=SIG_FRAMES_UPD{c}-intA;
        % SIG_FRAMES_UPD{c}=SIG_FRAMES_UPD{c};
        labelsUPD{c}=LABELS{1}(indexSIG);           % Clustering Labels
        NGroupsUPD{c}=length(unique(labelsUPD{c})); % N groups
        intA=intB;
        cummFrames=cummFrames+framesC;              % Cummulative Frames
    end
    NameCond=Condition_Names;
    disp('------------------>>Conditions Extracted.')
else
    NameCond=Names_Conditions;
end
%% Exporting to GEPHI and Ensemble Data per Condition
NG=0;
AuxC=0;
% Input: ColorState, from Clustering Analysis
for c=1:NCplot
    % Get raw data **************************************************
    R=Rasters{c}(:,IndexesActive);  % frames x cells
    ActNeu=find(sum(R,1)>0);        % Active Neurons in each Raster: ENS & NON ENS
    CAG=sum(R,2);                   % CAG in each Condition
    % Get labels 
    if UniRMutiE         % IF ALL EXPERIMETN in ONE Analysis
        labels=labelsUPD{c};            % Clustering Labels
        sigframes=SIG_FRAMES_UPD{c};    % Significative Frames
        NG=NGroupsUPD{c};               % N ensembles @ condition
        LENGHTRASTER{c}=size(R,1);
    else                 % Conditons were anayzed separately
        labels=LABELS{c};               % Clustering Labels
        sigframes=SIG_FRAMES{c};        % Significative Frames
        NG=NGroups{c};                  % N Ensembles @ condition
    end
    % Get Neural Ensembles Data **************************************
    if NG>0    %% IF THERE'RE ENSEMBLES in the RASTER
        % Re-Sorting in each Condition:
        [OrderOneCondition,Neurons_State_Cluster]=OrderClusters(labels,sigframes,R,NG);
        % CS=ColorState;
        % SORTING CONDITION
        Index_Ensemble=Indexes(OrderOneCondition);          % Neurons Label Raster
        % Hebbian Sequences 
        HebbSequence=Ensembles_Transitions(fs,labels,sigframes,ColorState,0); % ---> save
        % Save Network to Gephi**********************************
        ShNeuron = Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,ColorState,labels,Experiment);
        % ************************************************************
        % Features  of Ensembles ########################
        Ensembled_Labels{c}=labels;
        Ensemble_Threshold(c)=THR{c};
    else       %% Non-Ensembles
        ShNeuron=0;
        ActivityLength=0;
        % Features  of Ensembles ########################
        Ensembled_Labels{c}=[];
        Ensemble_Threshold(c)=THR{c};
    end
    % Update auxiliar variables
    if UniRMutiE
        AuxC=1;
    else
        AuxC=AuxC+NGroups{c};
    end

end
%% GET & SAVE NEURAL ENSEMBLE FEATURES
[Features_Ensemble,Features_Condition]=get_ensembles_features(R_Condition,Ensemble_Threshold,Ensembled_Labels,fs);

save_features_ensembles(Experiment,Names_Conditions,Features_Ensemble,Features_Condition);

%% Save SORTING,COORDINATES and COLOR ENSEMBLES
% Save this Only if it was Analyzed the whole enchilada
if CummFrames==TotalFrames
    checkname=1;
    while checkname==1
        ActualDir=pwd;
        SlashesIndx=find(ActualDir=='\');
        DefaultPath=[ActualDir(1:SlashesIndx(end)),'Processed Data'];
        if exist(DefaultPath,'dir')==0
            DefaultPath=pwd; % Current Diretory of Finder Spiker
        end
        [FileName,PathName] = uigetfile('*.mat',[' Pick the Analysis File ',Experiment],...
            'MultiSelect', 'off',DefaultPath);
        % dotindex=find(FileName=='.');
        if strcmp(FileName(1:end-4),Experiment)
            checkname=0;
            % SAVE DATA
            disp('>>Updating .mat File...')
            save([PathName,FileName],'XY_cluster','New_Order_Clustering',...
                'ColorState','-append');
            disp('>>Saving Analysis...')
            % Save Selected Analysis Variables
            for c=1:NC
                fprintf('>> Saving %s Analysis ... ',NameCond{c,1});
                save([PathName,FileName],VarNames{c},...
                '-append');
                fprintf('Done \n');
            end
            disp([Experiment,'   -> UPDATED (Ensembles Sorting Data)'])
        elseif FileName==0
            checkname=0;
            disp('ENSEMBLES SORTING DISCARDED')
            disp('NOT SAVED')
        else
            disp('Not the same Experiment!')
            disp('Try again!')
        end
    end
end

%% END OF THE WORLD