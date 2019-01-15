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
%  RASTER_Selected_Clean
%  XY_selected:         Coordinates Selected
%  fs
% Output
%  New_Order_Clustering:    Ensemble Sorting
%  XY_cluster:              Ensembled-Sorted Indexes
%  Experiment:              Concatenated Raster
%  Networks:                Gephi CSV files

%% Setup & Initialization
% addpath ..; % Functions from Signal Processing
[TotalCells,TotalFrames]=size(RASTER_Selected_Clean);
% NAMES PATCH
Experiment=Experiment(Experiment~='\');
NGroups={};
SIG_FRAMES={};
LABELS={};
Rasters={};
LENGHTRASTER={};
THR={};
Condition_Names={};
%% Select Number of COnditoins
NC = inputdlg('Enter Number of Analyzed Conditions:',...
         'Conditions', [1 50]);
NC= str2double(NC{:}); 

%% Read Result from NeuralNetworks Clustering Results Structure
for c=1:NC
    % Input Select Box
    WorkspaceVariables=who;
    % Filter Only with term _Analysis
    Filterby=['_Analysis'];
    [Nvar,~]=size(WorkspaceVariables);
    Variable_Index_Filtered=[];
    for nv=1:Nvar
        NamVar=WorkspaceVariables{nv};
        SameWords=ismember(Filterby,NamVar);
        if sum(SameWords)==length(Filterby)
            Variable_Index_Filtered=[Variable_Index_Filtered,nv];
        end
    end
    WorkspaceVariables=WorkspaceVariables(Variable_Index_Filtered);
    [index_var,index_CHECK] = listdlg('PromptString','Select a _Analysis Variable:',...
                    'SelectionMode','single',...
                    'ListString',WorkspaceVariables);
    % Get Structure Data Clustering Results
    if  index_CHECK>0
        VarNames{c}=WorkspaceVariables{index_var};
        NN_data=evalin('base',WorkspaceVariables{index_var}); % Data from NeuralNetworks
    end
    NamesFields=fieldnames(NN_data);
    if ismember('Clustering',NamesFields)   % Significant Coactivity
        THR{c,1}=NN_data.Peaks.Threshold;                        % Coactivity Threshold
        NGroups{c}=NN_data.Clustering.TotalStates;          % Number of Groups
        SIG_FRAMES{c}=find(NN_data.Peaks.Index);            % Significant COactivity Frames
        LABELS{c}=NN_data.Clustering.VectorStateIndex;      % Ensemble Labels
        % Ensemble_NBC_model{c}=NN_data.Classifier.Model;              % Naive Bayes Classifier
        % Ensemble_ValEerr{c}=NN_data.Classifier.ValidationError;       % Validation Error
    else                                    %  NO significant Coactivity
        NGroups{c}=0;      % Number of Groups
        SIG_FRAMES{c}=[];  % Significant COactivity Frames
        LABELS{c}=[];      % Ensemble Labels
    end
    % Get data from Strucutre
    Rasters{c}=NN_data.Data.Data;                       % Raster 
    LENGHTRASTER{c}=length(NN_data.Data.Data);          % Raster Length
    AuxName=WorkspaceVariables{index_var};              % Get Conditon: R_'Xcondition'
    indexes_name=find(AuxName=='_');
    Condition_Names{c}=AuxName(indexes_name(1)+1:indexes_name(end)-1);
end

%% Concatenation *****************************************************
% Initialize Variables
ExperimentRaster=[];    % Matrix Rx to concatenate Rasters
labels_frames=[];       % Ensemble Labels
signif_frames=[];       % Frames with Significant Coactivity
TotalNG=0;              % Total Nuber of Groups (N conditions->NG ensembles)
CummFrames=0;           % Cummulative Frames
CummGroups=0;           % Cummulative Groups
for i=1:NC
    % Raster
    ExperimentRaster=[ExperimentRaster;Rasters{i}]; % Frames x Cells
    % Significative Frames
    signif_frames=[makerowvector(signif_frames), makerowvector(SIG_FRAMES{i})+CummFrames];
    CummFrames=CummFrames+LENGHTRASTER{i};
    % Frame Labels
    labels_frames=[labels_frames;LABELS{i}+CummGroups];
    CummGroups=CummGroups+NGroups{i};
    % Group COunter:
    TotalNG=TotalNG+NGroups{i};
end
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

Indexes=1:length(XY_selectedClean);             % Indexes for AN
% IndexNeurons=Indexes;                       
Index_Ensemble=Indexes;
CoAc=sum(ExperimentRasterClean,2);              % Whole CoActivity
%% Sorting Clustering
[New_Order_Clustering,Neurons_State]=OrderClusters(labels_frames,signif_frames,ExperimentRasterClean,TotalNG);
% Indexes=sort(Indexes(New_Order_Clustering));    % useless
%% COLORMAP ENSEMBLES
ColorState=colormapensembles(TotalNG,NC,NGroups);
%% Plot Ensembles of Whole Raster ---------------------------------------------------------------------
%   Original ****************************************
OriginalExperiment=ExperimentRasterClean;
Plot_Raster_Ensembles(OriginalExperiment,Indexes,5,fs);                     % Disorted Raster
disp('Coloring Ensembles...')
Plot_State_Colors(labels_frames,signif_frames,ColorState,OriginalExperiment,fs,CoAc,Indexes);
disp('Coloring Ensembles Done.')
% plot_CAG_threshold(THR,R_Condition,fs)
plot_CAG_threshold(THR,LENGHTRASTER,fs)
if CummFrames==TotalFrames
    Label_Condition_Raster(Names_Conditions,R_Condition,fs);   % Labels
end
Figg=gcf; Figg.Name=['Neural Ensembles of ',Experiment];
%   Sorted *******************************************
Plot_Raster_Ensembles(OriginalExperiment,Indexes(New_Order_Clustering),1,fs);   % Sorted Raster
% Plot_State_Colors;
disp('Coloring Ensembles...')
Plot_State_Colors(labels_frames,signif_frames,ColorState,OriginalExperiment,fs,CoAc,Indexes(New_Order_Clustering));
disp('Coloring Ensembles Done.')
plot_CAG_threshold(THR,LENGHTRASTER,fs)
if CummFrames==TotalFrames
    Label_Condition_Raster(Names_Conditions,R_Condition,fs);   % Labels
end
Figg=gcf; Figg.Name=['Neural Ensembles (resorted) of ',Experiment];
% Ensemble Transitions HEBBIAN SEQUENCE **************
[ensemble_index_total]=Ensembles_Transitions(fs,labels_frames,signif_frames,ColorState,1);
XY_cluster=XY_selectedClean(New_Order_Clustering,:);    

%% Plot Ensembles in each Raster Condition (if NC>1)
% Initialize Output Features
HeadersDetailed={'Experiment','Condition','Ensemble','NeuronUse','TimeUse'};
AuxC=0;
NameSave={};
Ntransition=[];
RateSharedNeurons=[];
RateActivityLength=[];
Tdetailed=table();
EnsembleName={};
Ensembled_Raster={};
Ensembled_Labels={};
Ensemble_Threshold={};
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
%% Ploting and Exporting to GEPHI
% Checking if 1single raster has all conditions
NCplot=NC;
if CummFrames==TotalFrames && NC==1
        NCplot=length(Names_Conditions);
        disp('>>>> One Raster->All Conditions >>>>>')
        THR=repmat(THR,NCplot,1);
end
% EnsemblesCounter=0;
for c=1:NCplot
    NameSave{c,1}=Experiment(Experiment~='\');  % Experiment ID 
    R=Rasters{c}(:,IndexesActive);              % frames x cells
    if CummFrames==TotalFrames && NC==1 % IF ALL EXPERIMETN in ONE Analysis
        labels=labelsUPD{c};                    % Clustering Labels
        sigframes=SIG_FRAMES_UPD{c};            % Significative Frames
        % NG=numel(unique(labels_frames));
        NG=NGroupsUPD{c};                       % N ensembles @ condition
        LENGHTRASTER{c}=size(R,1);
    else                                % ELSE
        labels=LABELS{c};                   % Clustering Labels
        sigframes=SIG_FRAMES{c};            % Significative Frames
        NG=NGroups{c};                      % N Ensembles
    end
    ActNeu=find(sum(R,1)>0);   % Active Neurons in each Raster: ENS & NON ENS
    CoAc=sum(R,2);             % Coactivity Signal
    if NG>0     % IF THERE'RE ENSEMBLES
        % Re-Sorting in each Condition:
        [OrderOneCondition,Neurons_State_Cluster]=OrderClusters(labels,sigframes,R,NG);
        % ENSEMBLRE COLORS *******************************************
        if CummFrames==TotalFrames && NC==1
            % Included (deep) Purple:
            CS=ColorState;
        else
            CS = ColorState(AuxC+1:AuxC+NGroups{c},:);          
            CS=[CS;ColorState(end,:)];  % Plus (deep) Purple
        end
        % SORTING CONDITION
        Index_Ensemble=Indexes(OrderOneCondition);          % Neurons Label Raster
        % Plotting |
        Plot_Raster_Ensembles(R,Index_Ensemble,1,fs);       % Sorted Raster
        if CummFrames==TotalFrames && NC==1
            EnsembleFig=gcf; EnsembleFig.Name=Names_Conditions{c};
        else
            EnsembleFig=gcf; EnsembleFig.Name=Condition_Names{c};
        end
        disp('Coloring Ensembles...')
        Plot_State_Colors(labels,sigframes,CS,R,fs,CoAc,Index_Ensemble);
        disp('Coloring Ensembles Done.')
        plot_CAG_threshold(THR(c),LENGHTRASTER(c),fs)
        % check Name Condition: delete no alphanumeric symbosl
        okchar=[];
        for k=1:length(EnsembleFig.Name)
            if isalpha_num(EnsembleFig.Name(k))
                okchar=[okchar,k];
            end
        end
        % Create Ensemble Cells
        EnsembleName{c}=EnsembleFig.Name(okchar);
        Ensembled_Raster{c}=R;
        Ensembled_Labels{c}=labels;
        Ensemble_Threshold{c}=THR{c};
        % Hebbian Sequences 
        HebbSequence=Ensembles_Transitions(fs,labels,sigframes,CS,0); % ---> save
        % Cycles Reverberation Analysis: in waiting
        drawnow;
        % Save Network to Gephi**********************************
        % ShNeuron = NN2Gephi(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,ColorState,Experiment);
        % ShNeuron = Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,ColorState,labels,Experiment);
        ShNeuron = Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,CS,labels,Experiment);
        % ************************************************************
        ActivityLength=length(sigframes);
        EnsembleSorting{c}=OrderOneCondition;  % -------------->save
        Ntransition(c,1)=length(HebbSequence); % -------------->save
        % Ensembles Detailed Features
        NameCond={};
        EXP_ID={};
        NeuronUse=[];
        TimeUse=[];
        for i=1:NG
            EXP_ID{i,1}=Experiment(2:end);
            NeuronUse(i,1)=100*length(Neurons_State_Cluster{i})/length(ActNeu);
            TimeUse(i,1)=100*sum(labels==i)/length(sigframes);
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
            save([PathName,FileName],'EnsembleName','Ensembled_Raster',...
                'Ensembled_Labels','Ensemble_Threshold','UniRMutiE',...
                'ColorState','-append');
            % Ensemble_NBC_model,'Ensemble_ValEerr',...: not required
            disp('>>Save Ensemble Features...')
            disp('>>Saving Analysis...')
            % Save Analysis Variables
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

%% Nested Functions
% function plot
% end
%% END OF THE WORLD