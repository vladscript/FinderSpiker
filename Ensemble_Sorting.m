% Script to sort  and plot for Ensembles for Several Analized Conditions
% from Clustering Analysis from NeuralNetworks
% Input
%  NC:  How many Conditions: Input Dialogue
%  XY_Selected: Selected Active Coordinates
%  Raster_Analays Strucutre from NN: [ Frames x Cells ]
%  Strucutre from NN:   input from workspace:
%                       Raster
%                       Significative Frames
%                       Label of Ensembles
%  RASTER_Selected_Clean: Total SELECTED Frames
%  fs
% Output
%  New_Order_Clustering:    Ensemble Sorting
%  XY_cluster:              Ensembled-Sorted Indexes
%  Experiment:              Concatenated Raster
%  Networks:                Gephi CSV files

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

%% Plot Ensembles in each Raster Condition (if NC>1)
% Initialize Output Features
HeadersDetailed={'Experiment','Condition','Ensemble','NeuronUse','TimeUse'};
AuxC=0;
NameSave={};
Ntransition=[];
RateSharedNeurons=[];
RateActivityLength=[];
Tdetailed=table();
% EnsembleName={};
% Ensembled_Raster={};
% Ensembled_Labels={};
% Ensemble_Threshold={};
UniRMutiE=0;
%% IF ALL EXPERIMETN in ONE Analysis
if CummFrames==TotalFrames && NC==1
    disp('---------------------------------------------------------------')
    disp('>>>>>>>>>>>>>>>>   Whole Experiment   <<<<<<<<<<<<<<<<<<<<<<<<<')
    disp('---------------------------------------------------------------')
    UniRMutiE=1;
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
end
%% Exporting to GEPHI and Save Features per Condition
% Checking if 1single raster has all conditions
NCplot=NC;
if CummFrames==TotalFrames && NC==1
        NCplot=length(Names_Conditions);
        disp('>>>> One Raster->All Conditions >>>>>')
        THR=repmat(THR,NCplot,1);
end
% EnsemblesCounter=0;
for c=1:NCplot
    NameSave{c,1}=Experiment;           % Experiment ID 
    R=Rasters{c}(:,IndexesActive);      % frames x cells
    if CummFrames==TotalFrames && NC==1 % IF ALL EXPERIMETN in ONE Analysis
        labels=labelsUPD{c};            % Clustering Labels
        sigframes=SIG_FRAMES_UPD{c};    % Significative Frames
        % NG=numel(unique(labels_frames));
        NG=NGroupsUPD{c};                       % N ensembles @ condition
        LENGHTRASTER{c}=size(R,1);
    else                                % ELSE
        labels=LABELS{c};                   % Clustering Labels
        sigframes=SIG_FRAMES{c};            % Significative Frames
        NG=NGroups{c};                      % N Ensembles
    end
    ActNeu=find(sum(R,1)>0);   % Active Neurons in each Raster: ENS & NON ENS
    CAG=sum(R,2);             % Coactivity Signal
    if NG>0     % IF THERE'RE ENSEMBLES
        % Re-Sorting in each Condition:
        [OrderOneCondition,Neurons_State_Cluster]=OrderClusters(labels,sigframes,R,NG);
        % ENSEMBLRE COLORS *******************************************
        if CummFrames==TotalFrames && NC==1
            % Included (deep) Purple:
            CS=ColorState;
        else
            CS=ColorState(AuxC+1:AuxC+NGroups{c},:);          
            CS=[CS;ColorState(end,:)];  % Plus (deep) Purple
        end
        % SORTING CONDITION
        Index_Ensemble=Indexes(OrderOneCondition);          % Neurons Label Raster
        % Hebbian Sequences 
        HebbSequence=Ensembles_Transitions(fs,labels,sigframes,CS,0); % ---> save
        Ntransition(c,1)=numel(HebbSequence);
        % Cycles Reverberation Analysis: in waiting
        % drawnow;
        % Save Network to Gephi**********************************
        ShNeuron = Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,CS,labels,Experiment);
        % ************************************************************
        ActivityLength=length(sigframes);
        % EnsembleSorting{c}=OrderOneCondition;  % -------------->save
        % Ntransition(c,1)=length(HebbSequence); % -------------->save
        
        % Features  of Ensembles ########################
        NameCond={};
        EXP_ID={};
        NeuronUse=[];
        TimeUse=[];
        for i=1:NG
            EXP_ID{i,1}=Experiment;
            NeuronUse(i,1)=length(Neurons_State_Cluster{i})/length(ActNeu);
            TimeUse(i,1)=sum(labels==i)/length(sigframes);
            if CummFrames==TotalFrames && NC==1
                NameCond{i,1}=Names_Conditions{c};
                EnsembleList=(1:NG)';
            else
                NameCond{i,1}=Condition_Names{c};
                EnsembleList=(AuxC+1:AuxC+NG)';
            end
        end
        %Taux=table(repmat(Experiment(2:end),NG,1),NameCond,(1:NG)',NeuronUse,TimeUse);
        Taux=table(EXP_ID,NameCond,EnsembleList,NeuronUse,TimeUse);
        Taux.Properties.VariableNames=HeadersDetailed;
        Tdetailed=[Tdetailed;Taux];
    else            % Non-Ensembles
        ShNeuron=0;
        ActivityLength=0;
        Ntransition(c,1)=0;
        Plot_Raster_V(R,fs);
        if CummFrames==TotalFrames && NC==1
            EnsembleFig=gcf; EnsembleFig.Name=Names_Conditions{c};
        else
            EnsembleFig=gcf; EnsembleFig.Name=Condition_Names{c};
        end
        drawnow;
        % check Name Condition: delete no alphanumeric symbosl
        okchar=[];
        for k=1:length(EnsembleFig.Name)
            if isalpha_num(EnsembleFig.Name(k))
                okchar=[okchar,k];
            end
        end
        EnsembleName{c}=EnsembleFig.Name(okchar);
        Ensembled_Raster{c}=R;
        Ensembled_Labels{c}=[];
        Ensemble_Threshold{c}=THR{c};
    end
    % Ensembles Features
    RateSharedNeurons(c,1)=100*ShNeuron/length(ActNeu);
    if CummFrames==TotalFrames && NC==1
        RateActivityLength(c,1)=100*ActivityLength/LENGHTRASTER{1};
    else
        RateActivityLength(c,1)=100*ActivityLength/LENGHTRASTER{c};
    end

    % Update auxiliar variables
    if CummFrames==TotalFrames && NC==1
        AuxC=1;
    else
        AuxC=AuxC+NGroups{c};
    end

end

%% Save Ensembles Intel
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
                save([PathName,FileName],VarNames{c},...
                '-append');
            end
            disp([Experiment,'   -> UPDATED (Ensembles Data)'])
        elseif FileName==0
            checkname=0;
            disp('ENSEMBLES DISCARDED')
            disp('NOT SAVED')
        else
            disp('Not the same Experiment!')
            disp('Try again!')
        end
    end    
%% General Features Table *************************************************
HeadersFeatures={'Experiment','Condition','Ensembles','Transitions','Threshold','SharedNeurons'...
                    'ActivityDuration'};
if CummFrames==TotalFrames && NC==1
    Tensembles=table(NameSave,Names_Conditions,cell2mat(NGroupsUPD'),Ntransition,...
    THR,RateSharedNeurons,RateActivityLength,...
    'VariableNames',HeadersFeatures);
else
    Tensembles=table(NameSave,Condition_Names',cell2mat(NGroups'),Ntransition,...
    THR,RateSharedNeurons,RateActivityLength,...
    'VariableNames',HeadersFeatures);
end
disp(Tensembles);
disp(Tdetailed);

%% END OF THE WORLD