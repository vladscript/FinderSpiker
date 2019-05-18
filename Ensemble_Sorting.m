% Script to sort  and merge for Rasters for Several Analized Conditions
% from script get_bayes_ensembles and from GUI NeuralNetworks (by JP)
% 
% Requirement: load a .mat File from FinderSpiker
% 
% Input
%   Raster_Analysis      Data Structure of Ensemble Anlysis
%   From Worksapce:
%       RASTER_Selected_Clean
%       XY_Selected
%       fs
%       Experiment ID
% Output @ .mat File
%   New_Order_Clustering:    Sorting by Ensemble 
%   XY_cluster:              Ensembled-Sorted Indexes
%   Features_Ensemble
%   Features_Condition
%   ColorState
% Output @ .csv File
%   Neural General Ensemble Features:
%   Neural Detailed Ensemble Features:
%   Functional Network List of NODES & EDGES to import @ Gephi

%% Setup
Update_Directory;
%% Select Number of Analyzed Rasters ######
% NC = inputdlg('Enter Number of Analyzed Conditions:',...
%          'Conditions', [1 50]);
% NC= str2double(NC{:}); 
NC=numel(R_Condition);
TxtCnd{1,1}=['Experimental Conditions of ',num2str(Experiment),':'];
for n=1:NC
    TxtCnd{n+1,1}=Names_Conditions{n};
end
InitialMSG=msgbox(TxtCnd,'Experimental Conditions:');
% InitialMSG.Position=[450,340,125,130];
uiwait(InitialMSG);
delete(InitialMSG);

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

%% Select Sort of Color for Gephi Visualization
answer = questdlg('Make Node Colors for Gephi Visualization?', ...
	'Mix Color Ensembles', ...
	'MIX','EVEN','EVEN');


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
disp('>>Sorting Neurons')
if TotalNG<6
    [New_Order_Clustering,~]=OrderClusters(labels_frames,signif_frames,ExperimentRasterClean,TotalNG);
    re_sort=true;
else
    Nens=cell2mat(NGroups); n=1; oksort=true;
    re_sort=false;
    while oksort
        if Nens(n)<6
            oksort=false;
            fprintf('>>Re-sorting according to %s Ensembles\n',Names_Conditions{n});
            [New_Order_Clustering,~]=OrderClusters(labels_frames,signif_frames,ExperimentRasterClean,Nens(n));
            re_sort=true; 
        end
        n=n+1;
    end
    if ~re_sort
        New_Order_Clustering=Indexes;
        disp('>>No resort')
    end
end
disp('>>Neurons Sorted.')
% Coordinates
XY_cluster=XY_selectedClean(New_Order_Clustering,:); % Re-SORTED COORDINATES OF ENSEMBLES
% Indexes=sort(Indexes(New_Order_Clustering));    % useless
%% COLORMAP ENSEMBLES
% ColorState=colormapensembles(TotalNG,NC,NGroups);
ColorState=colormyensembles(NGroups);
% GUI to Choose Color SET
% Necessary: add CBREWER third-party function
% % % % % % % Nens=10;
% % % % % % % CT=cbrewer('qual','Set1',Nens);
% % % % % % % imagesc([1:Nens])
% % % % % % % colormap(CT)
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
if re_sort
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
end

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
NetworkCondition=cell(NCplot,1);
% Input: ColorState, from Clustering Analysis
for c=1:NCplot
    % Initialize Output
    % Max Links Between Neurons
    Network{c}.MaxSynLinks=0;
    % Max Links Between Neurons
    Network{c}.MaxCoupledPair=[0,0]; % Neurons A,B
    % Satitstics of Weights Connections (Synaptic Strength)
    Network{c}.SynStrengthStats=[0,0,0,0];
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
        switch answer
            case 'MIX'
                fNet=Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,ColorState,labels,Experiment,1);
            case 'EVEN'
                fNet=Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,ColorState,labels,Experiment);
            otherwise
                fNet=Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,ColorState,labels,Experiment);
        end
        %  Color Mode: Compensate:
        % fNet=Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,ColorState,labels,Experiment);
        %  Color Mode: MIXER:
        % fNet=Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,ColorState,labels,Experiment,1);
        
        
        NetworkCondition{c}=fNet;
        ShNeuron(c)=fNet.SharedNeruons/numel(ActNeu);
        % ************************************************************
        % Features  of Ensembles ########################
        Ensembled_Labels{c}=labels;
        Ensemble_Threshold(c)=THR{c};
    else       %% Non-Ensembles
        ShNeuron(c)=0;
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
Features_Condition.CoreNeurons=ShNeuron;
% Network Features
Features_Condition.Network=NetworkCondition;
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
            disp('>>Updating .mat with Neural Ensemble Analysis...')
            save([PathName,FileName],'XY_cluster','New_Order_Clustering',...
                                    'Features_Ensemble','Features_Condition',...
                                    'ColorState','-append');
            disp('>>Done.')
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
disp('>>END.')
%% END OF THE WORLD