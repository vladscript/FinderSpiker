Experiment=Experiment(Experiment~='\');     % NAMES PATCH
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
auxc=1;
for c=1:NC
    if ~isempty(index_var{c})
        % Get Structure Data Clustering Results
        VarNames{auxc}=WorkspaceVariables{index_var{c}}; % USED TO SAVE ANALYSIS VAR
        NN_data_Read{auxc}=evalin('base',WorkspaceVariables{index_var{c}}); % Data from WORKSPACE
        % Get Conditon Names from Variables
        AuxName=WorkspaceVariables{index_var{c}};              % Get Conditon: R_'Xcondition'
        indexes_name=find(AuxName=='_');
        Condition_Names{auxc}=AuxName(indexes_name(1)+1:indexes_name(end)-1);
        auxc=auxc+1;
    else
    end
end

%% Get DATA from Analysis Variables and CONCATENATE
% RASTERS ARE: [FRAMES x CELLS] -----------
[Rasters,LENGHTRASTER,NGroups,SIG_FRAMES,LABELS,THR,...
 ExperimentRaster,signif_frames,~,CummFrames,CummGroups,...
 TotalNG]=get_ensemble_intel(NN_data_Read);
clear NN_data_Read;
NC=numel(Rasters);

%% RE-SORT Ensembles labels
labels_frames=[];
offsetframes=0;
% signif_frames_INTER=[];
% labels_frames_INTER=[];
for c=1:NC
    sigfr=SIG_FRAMES{c}; % Starts @ 1 in each condition
    lbfr=LABELS{c};
    CAGc=sum(Rasters{c},2)';
    if ~isempty(lbfr)
        minEns=min(unique(lbfr))-1;
        [HB,HBtimes,HBinter]=Ensembles_Transitions(fs,lbfr-minEns,sigfr,CAGc,[],0,[],LENGHTRASTER{c});
        relbls=relabel_ensembles(lbfr-minEns,HB,'2-freq')+minEns;
    else
        minEns=0;
        HB=[]; HBtimes=[]; HBinter=[];
        relbls=[];
    end
    LABELS{c}=relbls;
    labels_frames=[labels_frames;relbls];
%     HebbLimits=HBinter+offsetframes;
%     HebbIndex=HB+minEns;
%     for hibb=1:numel(HB)
%         signif_frames_INTER=[signif_frames_INTER,HebbLimits(hibb,1):HebbLimits(hibb,2)];
%         labels_frames_INTER=[labels_frames_INTER,HebbIndex(hibb)*ones(1,HebbLimits(hibb,2)-HebbLimits(hibb,1)+1)];
%     end
%     offsetframes=LENGHTRASTER{c};
end

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
    while oksort && n<=numel(Nens)
        if Nens(n)<=6
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
% Original **************************************** *********************
OriginalExperiment=ExperimentRasterClean;
Plot_Raster_Ensembles(OriginalExperiment',fs,5,Indexes);  % Disorted Raster
disp('Coloring Ensembles...')
Plot_State_Colors(labels_frames,signif_frames,ColorState,OriginalExperiment,fs,CAG,Indexes);
disp('Coloring Ensembles Done.')
% plot_CAG_threshold(THR,R_Condition,fs)
plot_CAG_threshold(THR,LENGHTRASTER,fs)
% Labels
if CummFrames==TotalFrames
    Label_Condition_Raster(Names_Conditions,R_Condition,fs);   
else
    Label_Condition_Raster(Condition_Names,Rasters,fs);
end
Figg=gcf; Figg.Name=['Neural Ensembles of ',Experiment];
%% Ensemble Transitions HEBBIAN SEQUENCE ************** ALL-frame Details
Ensembles_Transitions(fs,labels_frames,signif_frames,CAG,ColorState,1,OriginalExperiment',LENGHTRASTER);
close(gcf); % Justo to color CAG  @ raster figure
% Sorted ******************************************************************
if re_sort
    Plot_Raster_Ensembles(OriginalExperiment',fs,1,Indexes(New_Order_Clustering));   % Sorted Raster
    % Plot_State_Colors;
    disp('Coloring Ensembles...')
    Plot_State_Colors(labels_frames,signif_frames,ColorState,OriginalExperiment,fs,CAG,Indexes(New_Order_Clustering));
    disp('Coloring Ensembles Done.')
    plot_CAG_threshold(THR,LENGHTRASTER,fs)
    if CummFrames==TotalFrames
        Label_Condition_Raster(Names_Conditions,R_Condition,fs);   % Labels
    else
        Label_Condition_Raster(Condition_Names,Rasters,fs);   % Labels
    end
    Figg=gcf; Figg.Name=['Neural Ensembles (resorted) of ',Experiment];
end

%% Ensemble Transitions HEBBIAN SEQUENCE ************** ALL-frame Details
Ensembles_Transitions(fs,labels_frames,signif_frames,CAG,ColorState,1,OriginalExperiment',LENGHTRASTER);