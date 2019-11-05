% Features of Neurons @ Ensembles
% It gets the Neurons that participate in each ensemble
% R(F)un after Clustering only
% INPUT:
%   R_Condition:        Cell of Rasters from Selected Raster
%   Ensemble_Threshold  Cell of vectors
%   Ensembled_Labels    Cell of vectors
%   fs                 sampling frequency
% OUTPUT
%   Ensembled_Neurons
%   Ensmeble_Features
function [Features_Ensemble,Features_Condition]=get_ensembles_features(R_Condition,Ensemble_Threshold,Ensembled_Labels,fs,LENGHTRASTER)
%% Setup
%SimMethod='hamming';
Load_Default_Clustering;
% Get Number Of Conditions
C=numel(R_Condition); % Number of Conditions
% Get Number of Different Ensmebles
Ensambles=[]; Nenscond=zeros(C,1);
for c=1:C
    Ensambles=[Ensambles;unique(Ensembled_Labels{c})];
    Nenscond(c)=numel(unique(Ensembled_Labels{c}));
end
Ensambles=unique(Ensambles);
Ne=numel(Ensambles); % TOTAL NUMBER OF ENSEMBLES (all experiment)
% Initialize OUTPUTS        Condition  x Ensembles
Ensembled_Neurons=cell(C,max(Nenscond));
Ensembles_Rate=zeros(C,max(Nenscond));
NeuronsOccupancy=zeros(C,max(Nenscond));
TimeOccupancy=zeros(C,max(Nenscond));
EnsCAGstats={};
% -> OUTPUT                 Condition
DunnIndex=zeros(1,C);
Transitions=cell(1,C);
Ntransitions=zeros(1,C);
Rate_Transitions=zeros(1,C);
Ratecycles=zeros(1,C);
CyclesTypes=zeros(3,C);
CAGstats=zeros(C,5);
% Thresholds=zeros(C,1);
ECV_Cond=zeros(C,1);
Model_Cond=cell(C,1);
MIV=zeros(C,1);
% TypeCycles(1)->Simple
% TypeCycles(2)->Closed
% TypeCycles(3)->Open
%% Main Loop to get Ensemble Features
for c=1:C
    %% DATA
    R=R_Condition{c};               % RASTER
    [AN,Frames]=size(R);            % Total Neurons [selected]
    if AN>Frames
        R=R';
        [AN,Frames]=size(R);            % Total Active Neurons [selected]
        fprintf('Transposed Matrix\n');
    end 
    RasterDuration=Frames/fs/60;    % MINUTES
    CAG=sum(R);                     % Co-Activity-Graphy
    %% CAG Statistics
    AUC=autocorr(CAG,1);                % CAG Autocorrelation Coefficient
    CAGstats(c,:)=AUC(2);
    Th=Ensemble_Threshold(c);               % CAG Threshold
    signif_frames=find(CAG>=Th);            % Significant Frames
    Ensembles_Labels=Ensembled_Labels{c};   % Labels each frame
    E=unique(Ensembled_Labels{c});          % Ensembles per condition
    %% Classification Error
    if ~isempty(Ensembles_Labels)
        [Model_Cond{c},ECV_Cond(c)]=Nbayes_Ensembles(R(:,signif_frames),Ensembles_Labels);
    else
        Model_Cond{c}={};
        ECV_Cond(c)=NaN;
    end
    %% HEBBIAN PATHWAY
    if isempty(signif_frames)
        HebbSequence=[];
        EnsembleTimes=[];
    else
        [HebbSequence,EnsembleTimes,EnsembleIntervals]=Ensembles_Transitions(fs,Ensembles_Labels,signif_frames,CAG,[],0,LENGHTRASTER{c});
    end
    %% Ensembles DUration & Inter Ensemble Interval
    % Inter-Eensemble-Interval & Ensemble Duration*********************
    % FROM HEBBIAN PATH 
    EDs=(diff(EnsembleIntervals')+1)/fs; % [SAMPLES]
    if size(EnsembleIntervals,1)>1
%         for i=2:size(EnsembleIntervals,1)
%             IEIs(1,i-1)=(EnsembleIntervals(i,1)-EnsembleIntervals(i-1,2))/fs;
%         end
        IEIs=diff(EnsembleTimes)/fs;
    else
        IEIs=[];
    end
    % INTER (ALL) ENSEMBLES INTERVAL
    allIEIsExp{c}=IEIs;   % [SECONDS]
    % (ALL) ENSEMBLES DURATION
    allEDsExp{c}=EDs;     % [SECONDS]
    %% EACH ENSEMBLE FEATURES
    MaxIntraVec=zeros(1,numel(E));
    for e=1:numel(E)
        fprintf('>> Condition %i, Ensemble %i of %i \n',c,e,numel(E));
        frames_ensemble=signif_frames(Ensembles_Labels==E(e));
        TimeOccupancy(c,e)=numel(frames_ensemble)/numel(signif_frames);
        EnsembleActivations=numel(HebbSequence(HebbSequence==E(e)));
        % Ouput Measure Features by ENSEMBLE
        Ensembles_Rate(c,e)=EnsembleActivations/RasterDuration; %[act/min]
        Ensembled_Neurons{c,e}=find(sum(R(:,frames_ensemble),2));
        % Neurons@Ensemble / Active Neurons @ Condition
        NeuronsOccupancy(c,e)=numel(Ensembled_Neurons{c,e})/sum(sum(R,2)>0); 
        % MORE FEATURES
        Rcluster=R(:,frames_ensemble); % Cells x Frames
        CAGcluster=sum(Rcluster);
        if numel(CAGcluster)>1
            AUC=autocorr(CAGcluster,1); % Autocorrelation Coeffcient
            EnsCAGstats{c,e}=[AUC(2),mean(CAGcluster),mode(CAGcluster),median(CAGcluster),var(CAGcluster),skewness(CAGcluster),kurtosis(CAGcluster)];
            MaxIntraVec(c)=max(pdist(Rcluster',SimMethod));
        else
            EnsCAGstats{c,e}=[0,mean(CAGcluster),mode(CAGcluster),median(CAGcluster),var(CAGcluster),skewness(CAGcluster),kurtosis(CAGcluster)];
            MaxIntraVec(c)=0;
        
        end
        eEDs=EDs(HebbSequence==E(e)); % [seconds]
        eIEIs=diff(EnsembleTimes(HebbSequence==E(e)))/fs; % [samples/fs=seconds]
        % Interval (same) Ensemble Inter
        IEIstats{c,e}=[mean(eIEIs),mode(eIEIs),median(eIEIs),var(eIEIs),skewness(eIEIs),kurtosis(eIEIs)];
        % Ensemble Duration
        EDstats{c,e}=[mean(eEDs),mode(eEDs),median(eEDs),var(eEDs),skewness(eEDs),kurtosis(eEDs)];
        IEIsExp{c}=eIEIs;   % [SECONDS]
        EDsExp{c}=eEDs;     % [SECONDS]
        
    end
    if isempty(MaxIntraVec)
        MIV(c)=NaN;
    else
        MIV(c)=max(MaxIntraVec);
    end
    %% ENSEMBLES SETs FEATURES
    % CONDITION FEATURES ********************************
    % Dunn's Index (sort of):
    % How Separater Cluster divided how Divided Intra Clusters
    % <1 More Distance Intra Vectors than Intra Clusters-> Bad Clustering
    % >1 More Distance Intra Clusters than Intra Vectors-> Good Clustering
    NeuroClusters=zeros(AN,numel(E));
    for e=1:length(E)
        NeuroClusters(Ensembled_Neurons{c,e},e)=1;
    end
    Dhamm=pdist(NeuroClusters',SimMethod); % percentage of different neurons
    if isempty(Dhamm); Dhamm=0; end;
    if max(MaxIntraVec)>0
        DunnIndex(c)=min(Dhamm)/max(MaxIntraVec); % min distance among ensembles divided by maximum length of ensembles
    else
        DunnIndex(c)=0; % min distance among ensembles divided by maximum length of ensembles
    end
    
    %% TRANSITIONS AND CYCLES
    ET=HebbSequence;
    Transitions{c}=ET;
    Ntransitions(c)=numel(ET)-1;
    Rate_Transitions(c)=(numel(ET)-1)/RasterDuration; % Transitions per MINUTE
    % Alternate or Reactivate Proportion of the Sequence
    Reactivations=0;
    Alternations=0;
    AltIndex(c)=0;
    ReaIndex(c)=0;
    if numel(HebbSequence)>1
        for e=1:numel(HebbSequence)-1
             if HebbSequence(e)==HebbSequence(e+1)
                 Reactivations=Reactivations+1;
             else
                 Alternations=Alternations+1;
             end
        end
        AltIndex(c)=Alternations/(numel(HebbSequence)-1);
        ReaIndex(c)=Reactivations/(numel(HebbSequence)-1);
    end
    
    % Euler Cycles of Ensembles [REVERBERATION] return 
    % to any given ensemble 
    % (after activate all ensembles)
    [TableCycles,TypeCycles]=eulercycles_hebbseq(ET,E);
    CyclesTypes(:,c)=TypeCycles;
    Ratecycles(c)=sum(TypeCycles)/RasterDuration; % [cycles/min]
    % Timing of the Cycles #########################################
    Ncy=size(TableCycles,1);
    Starkts=[]; % Saver Start Index of THe Cycles
    if Ncy>0
        for n=1:Ncy
            Cycle=TableCycles{n,2};
            searchin=true;
            nn=1;
            while searchin
                SegmentDiffE=ET(nn:nn-1+numel(Cycle))-Cycle;
                if sum(SegmentDiffE==0)==numel(Cycle)&&~ismember(nn,Starkts)
                    TableCycles{n,3}=EnsembleTimes(nn)/fs; % [SECONDS]
                    TableCycles{n,4}=EnsembleTimes(nn-1+numel(Cycle))/fs; % [SECONDS]
                    Starkts=[Starkts,nn];
                    disp('>>Found Cycle')
                    %ETcopy(nn)=-666; % unmark!
                    searchin=false;
                end
                nn=nn+1;
            end
        end
    end
    TableCyclesCell{c}=TableCycles;
    %     TableCYcles:
    % Type      Sequence     ConditioinStart[s] ConditionEnd[s]:
    % 'simple'  1,2,3,1      10                     12    
    
    if c<C; disp('Next Condition'); end;
end
%% Cross Simmilar Neural Ensembles Are 
NeuralClusters=zeros(AN,Ne);
col=1;
for ncond=1:C
    for nens=1:Nenscond(ncond)
        NeuralClusters(Ensembled_Neurons{ncond,nens},col)=1;
        col=col+1;
    end
end
% CrossEnsmebleSimm=ones(Ne);
CrossEnsmebleSimm=1-squareform( pdist(NeuralClusters',SimMethod) );


%% OUTPUTS ***************************************************************

%                       Vectors of Condition x 1
% ABOUT the CLUSTERING
Features_Condition.Nenscond=Nenscond;
Features_Condition.Threshold=Ensemble_Threshold;
Features_Condition.Dunn=DunnIndex;
Features_Condition.MIV=MIV;
Features_Condition.ECV_Cond=ECV_Cond;
% HEBBIAN SEQUENCE
Features_Condition.HebbianSeq=Transitions;
Features_Condition.HebbianDurations=allEDsExp; % *
Features_Condition.HebbianInterval=allIEIsExp; % *
% ABOUT the ALTERNANCE & REVERBERATION
Features_Condition.RateTrans=Rate_Transitions;
Features_Condition.RateCycles=Ratecycles;
Features_Condition.CyclesType=CyclesTypes;
Features_Condition.AltIndx=AltIndex;
Features_Condition.ReaIndx=ReaIndex;
% ABOUT CO-ACTIVITY GRAPHY
Features_Condition.CAGstats=CAGstats;

%                                   MAT File:
Features_Condition.Model_Cond=Model_Cond;
Features_Condition.CrossEnsmebleSimm=CrossEnsmebleSimm;
Features_Condition.CyclesTable=TableCyclesCell; % Type and Sequence

% Domminance:
% Dominance_Ensemble=NeuronsOccupancy.*Ensembles_Rate;
Dominance_Ensemble=NeuronsOccupancy.*TimeOccupancy;
%                           Matrices of Conditon x Ensembles
Features_Ensemble.Neurons=Ensembled_Neurons;
% Time and Size of Ensembles
Features_Ensemble.NeuronsOccupancy=NeuronsOccupancy;
Features_Ensemble.TimeOccupancy=TimeOccupancy;
Features_Ensemble.Dominance=Dominance_Ensemble;
Features_Ensemble.Rate=Ensembles_Rate;
% Ensemble CAG features
Features_Ensemble.EnsCAGstats=EnsCAGstats;
% Ensemble Time features
Features_Ensemble.IEIstats=IEIstats; % [seconds]
Features_Ensemble.EDstats=EDstats;   % [seconds]

%                                   MAT File
Features_Ensemble.IEIsExp=IEIsExp;   % [seconds]
Features_Ensemble.EDsExp=EDsExp;     % [seconds]

disp('>>Feature Extraction of Neural Ensemble: Done.')
%% disp('END OF THE WORLD')