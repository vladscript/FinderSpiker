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
function [Features_Ensemble,Features_Condition]=get_ensembles_features(R_Condition,Ensemble_Threshold,Ensembled_Labels,fs)
%% Setup
SimMethod='hamming';
% Get Number Of Conditions
C=numel(R_Condition); % Number of Conditions
% Get Number of Different Ensmebles
Ensambles=[];
for c=1:C
    Ensambles=[Ensambles;unique(Ensembled_Labels{c})];
end
Ensambles=unique(Ensambles);
Ne=numel(Ensambles); % TOTAL NUMBER OF ENSEMBLES (all experiment)

% Initialize Output: Condition  x Ensembles
Ensembled_Neurons=cell(C,Ne);   % -> OUTPUT per ENSEMBLE
Ensembles_Rate=zeros(C,Ne);     % -> OUTPUT per ENSEMBLE
NeuronsOccupancy=zeros(C,Ne);   % -> OUTPUT per ENSEMBLE
% MORE FEATURES ##########################################
DunnIndex=zeros(1,C);                % -> OUTPUT per Condition
HebbianSequence=cell(1,C);      
Transitions=cell(1,C);
Ntransitions=zeros(1,C);        
Rate_Transitions=zeros(1,C);    % -> OUTPUT per Condition
Ratecycles=zeros(1,C);          % -> OUTPUT per Condition
CyclesTypes=zeros(3,C);         % -> OUTPUT per Condition
% TypeCycles(1)->Simple
% TypeCycles(2)->Closed
% TypeCycles(3)->Open
% Ensemble_Features=cell(C,Ne);

%% Main Loop

for c=1:C
    %% DATA
    R=R_Condition{c};               % RASTER
    [AN,Frames]=size(R);            % Total Active Neurons [selected]
    RasterDuration=Frames/fs/60;    % MINUTES
    CAG=sum(R);                     % Co-Activity-Graphy
    % CAG Statistics
    AUC=autocorr(CAG,1); % Autocorrelation Coeffcient
    CAGstats=[AUC(2),mean(CAG),var(CAG),skewness(CAG),kurtosis(CAG)];
    Th=Ensemble_Threshold{c};       % CAG Threshold
    signif_frames=find(CAG>=Th);    % Significatn Frames
    Ensembles_Labels=Ensembled_Labels{c}; % Labels each frame
    if numel(signif_frames)==numel(Ensembles_Labels)
        E=unique(Ensembled_Labels{c}); % Ensambles per condition
        MaxIntraVec=zeros(1,numel(E));
        EnsCAGstats=zeros(numel(E),5);
        %% EACH ENSEMBLE
        for e=1:numel(E)
            fprintf('>> Condition %i, Ensemble %i of %i \n',c,e,numel(E));
            frames_ensemble=signif_frames(Ensembles_Labels==E(e));
            TimeOccupancy=numel(frames_ensemble)/numel(signif_frames);
            EnsembleActivations=numel(find(diff(frames_ensemble)>1));
            % Ouput Measure Features by ENSEMBLE
            Ensembles_Rate(c,e)=EnsembleActivations/RasterDuration;
            Ensembled_Neurons{c,e}=find(sum(R(:,frames_ensemble),2));
            % Output Indexes
            NeuronsOccupancy(c,e)=numel(Ensembled_Neurons{c,e})/AN;
            % MORE FEATURES
            Rcluster=R(:,frames_ensemble); % Cells x Frames
            CAGcluster=sum(Rcluster);
            AUC=autocorr(CAGcluster,1); % Autocorrelation Coeffcient
            EnsCAGstats(e,:)=[AUC(2),mean(CAGcluster),var(CAGcluster),skewness(CAGcluster),kurtosis(CAGcluster)];
            MaxIntraVec(e)=max(pdist(Rcluster',SimMethod));
            % Inter-Eensemble-Interval & Ensemble Duration
            IEIall=diff(frames_ensemble);
            IEI=IEIall(IEIall>1);
            disp('butt')
        end
    else
        disp('>error<') % likely impossible this ever happens
    end
    %% ENSEMBLES SETs FEATURES
    % [CAGstats,EnsCAGstats]-> SAVE as row in the csv
    NeuroClusters=zeros(AN,numel(E));
    for e=1:length(E)
        NeuroClusters(Ensembled_Neurons{c,e},e)=1;
    end
    % CONDITION FEATURES ********************************
    % Dunn's Index (sort of):
    % How Separater Cluster divided how Divided Intra Clusters
    % <1 More Distance Intra Vectors than Intra Clusters-> Bad Clustering
    % >1 More Distance Intra Clusters than Intra Vectors-> Good Clustering
    Dhamm=pdist(NeuroClusters',SimMethod); % percentage of different neurons
    if isempty(Dhamm); Dhamm=0; end;
    if max(MaxIntraVec)>0
        DunnIndex(c)=min(Dhamm)/max(MaxIntraVec); % min distance among ensembles divided by maximum length of ensembles
    else
        DunnIndex(c)=0; % min distance among ensembles divided by maximum length of ensembles
    end
    % MORE FEATURES #########
    
    % Hebbian Sequence
    HS=Ensembles_Labels(diff(signif_frames)>1);
    HebbianSequence{c}=HS;
    % Transitions: Ensmbles Change deactivate and activate [ALTERNANCE]
    ET= HS(diff([HS;0])~=0);
    Transitions{c}=ET;
    Ntransitions(c)=numel(ET);
    Rate_Transitions(c)=numel(ET)/RasterDuration; % Transitions per MINUTE
    % Cycles of Ensembles [REVERBERATION] return to a given ensemble (after activate all ensembles)
    TypeCycles=zeros(3,1);
    tcounter=[]; t=1;
    Tremaining=1;
    while and(~isempty(Tremaining),~isempty(ET))
        ActualEnsemble=ET(t);
        Cy=find(ET(t+1:end)==ActualEnsemble);
        if ~isempty(Cy)
            auxt=t;
            i=1;
            while i<=length(Cy)
                Cycle=ET(auxt:t+Cy(i));
                % Check what kind of Hebbian Path: only for Cycle with all
                % Active Ensemables
                % Simple
                if and(numel(Cycle)==length(E)+1,numel(unique(Cycle))==numel(E))
                    disp(Cycle')
                    disp('Simple')
                    TypeCycles(1)=TypeCycles(1)+1;
                    tcounter=[tcounter;t;Cy(1:i)+t];
                    auxt=t+Cy(i);
                    i=i+1;
                elseif numel(unique(Cycle))==numel(E)
                    disp(Cycle')
                    CycleMat=zeros(max(E));
                    for j=1:length(Cycle)-1
                        CycleMat(Cycle(j),Cycle(j+1))=CycleMat(Cycle(j),Cycle(j+1))+1;
                        % ASS(j,:)=[Cycle(j),Cycle(j+1)]
                    end
                    % Check for upper and lower triangles
                    CS=sum(triu(CycleMat)+triu(CycleMat'));
                    % If there were sequences from 2 to end
                    if sum(CS(2:end)>0)==numel(E)-1
                        disp('Closed')
                        TypeCycles(2)=TypeCycles(2)+1;
                        tcounter=[tcounter;t;Cy(1:i)+t];
                        auxt=t+Cy(i);
                        i=i+1;
                    else
                        disp('Open')
                        TypeCycles(3)=TypeCycles(3)+1;
                        tcounter=[tcounter;t;Cy(1:i)+t];
                        auxt=t+Cy(i);
                        i=i+1;
                    end
                else
                    tcounter=[tcounter;t;Cy(1:i)+t];
                    i=i+1;
                    auxt=t;
                end
                tcounter=unique(tcounter);
            end
        else
            Cycle=[];
            tcounter=unique([tcounter;t]);
        end
        Tremaining=setdiff(1:numel(ET)-1,tcounter);
        if ~isempty(Tremaining)
            t=Tremaining(1);
        end
    end
    CyclesTypes(:,c)=TypeCycles;
    Ratecycles(c)=sum(TypeCycles)/RasterDuration;
    if c<C; disp('Next Condition'); end;
end
%% Rate of Change of Neurons in Ensembles
% ONly for C>1 conditions
if C>1
    RateNeuronsChanges=zeros(C-1,Ne);
    for c=2:C
        for n=1:Ne
            ActualVector=zeros(AN,1);
            PreVector=zeros(AN,1);
            PreVector(Ensembled_Neurons{c-1,n})=1;
            ActualVector(Ensembled_Neurons{c,n})=1;
            PrevailedVector=PreVector.*ActualVector;
            if sum(ActualVector)<sum(PreVector)
                RateNeuronsChanges(c-1,n)=sum(PrevailedVector)/sum(PreVector);
                disp('Reduced Ensemble')
            else
                RateNeuronsChanges(c-1,n)=(sum(ActualVector)-sum(PreVector))/sum(PrevailedVector);
                if RateNeuronsChanges(c-1,n)<1
                    RateNeuronsChanges(c-1,n)=RateNeuronsChanges(c-1,n)+1;
                end
                disp('Enhanced Ensemble')
            end
        end
    end
else
    RateNeuronsChanges=zeros(1,Ne);
end
Dominance_Ensemble=NeuronsOccupancy.*Ensembles_Rate;
[~,Ensemble_Dominant]=max(Dominance_Ensemble')
%% OUTPUT **********************************************************
Features_Ensemble.Neurons=Ensembled_Neurons;
Features_Ensemble.Rate=Ensembles_Rate;
Features_Ensemble.Dominance=Dominance_Ensemble;
% Features_Ensemble.Threshold=Dominance_Ensemble;
% Features_Ensemble.AllEnsemble=Dominance_Ensemble;
Features_Condition.Dunn=DunnIndex;
Features_Condition.RateTrans=Rate_Transitions;
Features_Condition.RateCycles=Ratecycles;
Features_Condition.CyclesType=CyclesTypes;
% Features_Condition.Ensemble_Dominant=Ensemble_Dominant;
% disp('END OF THE WORLD')

