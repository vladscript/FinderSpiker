%% Setup
Update_Directory;
NC=numel(R_Condition);
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
%%
answer='YES';

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
ColorState=colormapensembles(TotalNG,NC,NGroups);

%% Plot Ensembles of Whole Raster ---------------------------------------------------------------------
Experiment=Experiment(Experiment~='\');     % NAMES PATCH
%   Original ****************************************
OriginalExperiment=ExperimentRasterClean;

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
            case 'YES'
                fNet=Get_Gephi_Network(R,XY_selectedClean,Neurons_State_Cluster,OrderOneCondition,ColorState,labels,Experiment,1);
            case 'NO'
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