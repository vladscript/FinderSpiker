%% Save DATA: Neural Ensembles Features
% .mat file:
%   Features_Ensemble
%   Features_Condition
% .CSV file
%   N-Ensembles,Dunns Index,RateOfTransitions,RateOfCycles,DominantEnsemble
%   N SimpleCycles, N-ClosedCycles , N-OpenedCycles
%   Ensemble i Rate,...
%   Ensemble i Dominance=%NeuronsOccupance * Rate,...
function save_features_ensembles(Experiment,Names_Conditions,Features_Ensemble,Features_Condition)
% Setup
% Saving Directory: one above where Finder Spiker is..
Experiment=Experiment(Experiment~='\');     % NAMES PATCH
FileDirSave=pwd;
slashes=find(FileDirSave=='\');
FileDirSave=FileDirSave(1:slashes(end));

%% SAVE OUTPUT DATASET (.m file)
% checkname=1;
% while checkname==1
%     DefaultPath=[FileDirSave,'Processed Data'];
%     if exist(DefaultPath,'dir')==0
%         DefaultPath=pwd; % Current Diretory of MATLAB
%     end
%     [FileName,PathName] = uigetfile('*.mat',[' Pick the Analysis File ',Experiment],...
%         'MultiSelect', 'off',DefaultPath);
%     dotindex=find(FileName=='.');
%     if strcmp(FileName(1:end-4),Experiment)
%         checkname=0;
%         % SAVE DATA
%         disp('>>Updating .mat file...')
%         save([PathName,FileName],'Features_Ensemble','Features_Condition',...
%             '-append');
%         disp([Experiment,'   -> UPDATED (Ensembles Features)'])
%     elseif FileName==0
%         checkname=0;
%         disp('*************DISCARDED************')
%     else
%         disp('Not the same Experiment!')
%         disp('Try again!')
%     end
% end    
%% SAVE CSV FILES
% Direcotry Name
NameDir='Ensemble Features\';
% Number of Condition
if iscell(Names_Conditions)
    C=numel(Names_Conditions);
else
    C=1;
end
%% Save Files for each COndition
HeadersFeaturesCondition={'Nensembles','Threshold','Dunns','MaxIntraVec','ClassError',...
                     'RateTrans','RateCycles','SimpleCyclesP','ClosedCyclesP','OpenedCyclesP',...
                     'CAGauc','CoreSize','MaxSynLinks','MaxConn_A','MaxConn_B',...
                     'EDMean','EDMode','EDMedian','EDVar','EDSkew','EDKurt',...
                     'IEIMean','IEIMode','IEIMedian','IEIVar','IEISkew','IEIKurt',...
                     'SynWeigthMean','SynWeigthMode','SynWeightMedian','SynWeigthVar','SynWeigthSkew','SynWeigthKurt',...
                     'AlternativeIndx','RecurrentIndx'};
HeadersFeaturesEnsembles={'NeuronsRation','Time','Dominance','Rate',...
     'ensCAGauc','ensCAGmean','ensCAGmode','ensCAGmedian','ensCAGvar','ensCAGskew','ensCAGkurt'...
     'IEImean','IEImode','IEImedian','IEIvar','IEIskew','IEIkurt'...
     'EDmean','EDmode','EDmedian','EDvar','EDskew','EDkurt'};
% For each Ensembles thera rows are : Ens1; Ens2; ... Ensj

for c=1:C
    fprintf('>> Creating Table Ensemble Feautures of %s \n',Names_Conditions{c});
    %% General INTEL about Neurla Ensembles
    % Features Columns of the Table ***************************************
    Name=Names_Conditions{c};
    NE=Features_Condition.Nenscond(c);                  % N ensembles
    Th=Features_Condition.Threshold(c);                 % Threshold
    DunnIndx=Features_Condition.Dunn(c);                % Dunns Index
    MaxIntraVec=Features_Condition.MIV(c);              % Max Distance Intra Vectors
    ClassError=Features_Condition.ECV_Cond(c);          % Classification Error
    RateTran=Features_Condition.RateTrans(c);           % Rate Transitions
    RateCycl=Features_Condition.RateCycles(c);          % Rate Cycles
    EDs=Features_Condition.HebbianDurations{c};         % Ensmeble Durations
    IEIs=Features_Condition.HebbianInterval{c};         % Ensmeble Intervals
    % [%] Simple Cycles
    Portion_simple=100*Features_Condition.CyclesType(1,c)/sum(Features_Condition.CyclesType(:,c));
    if isnan(Portion_simple); Portion_simple=0; end;
    % [%] Closed Cycles
    Portion_closed=100*Features_Condition.CyclesType(2,c)/sum(Features_Condition.CyclesType(:,c));
    if isnan(Portion_closed); Portion_closed=0; end;
    % [%] Opened Cycles
    Portion_opened=100*Features_Condition.CyclesType(3,c)/sum(Features_Condition.CyclesType(:,c));
    if isnan(Portion_opened); Portion_opened=0; end;
    CAGauc=Features_Condition.CAGstats(c,1);            % CAG AutoCorrelation Coefficient
    CoreSize=Features_Condition.CoreNeurons(c);         % Ratio of neurons in all Ensembles
    % Max Links Between Neurons
    if ~isempty(Features_Condition.Network{c})
        MaxSynLinks=Features_Condition.Network{c}.MaxSynLinks;
    else
        Features_Condition.Network{c}.MaxSynLinks=0;
        MaxSynLinks=Features_Condition.Network{c}.MaxSynLinks;
        Features_Condition.Network{c}.SynStrengthStats(1)=0;  % mean
        Features_Condition.Network{c}.SynStrengthStats(2)=0;  % mode
        Features_Condition.Network{c}.SynStrengthStats(3)=0;  % median
        Features_Condition.Network{c}.SynStrengthStats(4)=0;  % variance
        Features_Condition.Network{c}.SynStrengthStats(5)=0;  % skewness
        Features_Condition.Network{c}.SynStrengthStats(6)=0;  % kurtosis
    end
    % Max Links Between Neurons
    if MaxSynLinks>0
        MaxConnA=Features_Condition.Network{c}.MaxCoupledPair(1); % Neuron A
        MaxConnB=Features_Condition.Network{c}.MaxCoupledPair(2); % Neuron B
    else
        MaxConnA=0; % Neuron A
        MaxConnB=0; % Neuron B
    end
    % Satitstics of Weights Connections (Synaptic Strength)
%     from get_gephi[...] :
%       fNet.SynStrengthStats=[mean(WEIGHT),mode(WEIGHT),median(WEIGHT),var(WEIGHT),skewness(WEIGHT),kurtosis(WEIGHT)];
    SynMean=Features_Condition.Network{c}.SynStrengthStats(1);  % mean
    SynMode=Features_Condition.Network{c}.SynStrengthStats(2);  % mode
    SynMedian=Features_Condition.Network{c}.SynStrengthStats(3);% median
    SynVar=Features_Condition.Network{c}.SynStrengthStats(4);   % variance
    SynSkew=Features_Condition.Network{c}.SynStrengthStats(5);  % skewness
    SynKurt=Features_Condition.Network{c}.SynStrengthStats(6);  % kurtosis
    % Alternative Index
    AlternationIndex=Features_Condition.AltIndx(c);
    % Recurrent Index
    RecurrentIndex=Features_Condition.ReaIndx(c);
    % ********************************************************************
    % Create Table
    Tensemblesfeatures=table(NE,Th,DunnIndx,MaxIntraVec,ClassError,...
        RateTran,RateCycl,Portion_simple,Portion_closed,Portion_opened,...
        CAGauc,CoreSize,MaxSynLinks,MaxConnA,MaxConnB,...
        mean(EDs),mode(EDs),median(EDs),var(EDs),skewness(EDs),kurtosis(EDs),...
        mean(IEIs),mode(IEIs),median(IEIs),var(IEIs),skewness(IEIs),kurtosis(IEIs),...
        SynMean,SynMode,SynMedian,SynVar,SynSkew,SynKurt,...
        AlternationIndex,RecurrentIndex);
    Tensemblesfeatures.Properties.VariableNames=HeadersFeaturesCondition;
    % Save CSV
    if ~isdir([FileDirSave,NameDir])
        mkdir([FileDirSave,NameDir]);
        disp('>>Directory >Ensemble Features< created')
    end
    disp('>>Saving...') 
    writetable(Tensemblesfeatures,[FileDirSave,NameDir,Experiment,'_',Name,'.csv'],...
        'Delimiter',',','QuoteStrings',true);
    disp(['>>Saved at /Ensemble Features: ',Experiment,'-',Names_Conditions{c}])
    %% Detailed Ensembles INTEL ###########################################
    % Features For each Ensemble (rows of the table)
    % [%] Active Neurons (observed) [column vectors!]
    NeuronsOccupancy=Features_Ensemble.NeuronsOccupancy(c,1:NE)'; 
    % [%] Active Time (observed)
    TimeOccupancy=Features_Ensemble.TimeOccupancy(c,1:NE)';
    % Dominance Index
    DominanceIndex=Features_Ensemble.Dominance(c,1:NE)';
    % Ensemble Rate [act/min]
    EnseRate=Features_Ensemble.Rate(c,1:NE)';
    % *Initialize* sizes right
    EnsCAGauc=zeros(NE,1);EnsCAGmean=EnsCAGauc; EnsCAGmode=EnsCAGauc;
    EnsCAGmedian=EnsCAGauc; EnsCAGvar=EnsCAGauc;
    EnsCAGskew=EnsCAGauc; EnsCAGkurt=EnsCAGauc; 
    EnsIEImean=EnsCAGauc; EnsIEImode=EnsCAGauc; EnsIEImedian=EnsCAGauc;
    EnsIEIvar=EnsCAGauc; EnsIEIskew=EnsCAGauc; EnsIEIkurt=EnsCAGauc;
    EnsEDmean=EnsCAGauc; EnsEDmode=EnsCAGauc; EnsEDmedian=EnsCAGauc;
    EnsEDvar=EnsCAGauc; EnsEDskew=EnsCAGauc;
    EnsEDkurt=EnsCAGauc;
    for n=1:NE
        % Ensemble Co-Activity-Graphy ACC *******************
        EnsCAGauc(n,1)=Features_Ensemble.EnsCAGstats{c,n}(1);
        % Ensemble Co-Activity-Graphy Mean
        EnsCAGmean(n,1)=Features_Ensemble.EnsCAGstats{c,n}(2);
        % Ensemble Co-Activity-Graphy Mode
        EnsCAGmode(n,1)=Features_Ensemble.EnsCAGstats{c,n}(3);
        % Ensemble Co-Activity-Graphy Median
        EnsCAGmedian(n,1)=Features_Ensemble.EnsCAGstats{c,n}(4);
        % Ensemble Co-Activity-Graphy Variance
        EnsCAGvar(n,1)=Features_Ensemble.EnsCAGstats{c,n}(5);
        % Ensemble Co-Activity-Graphy Skewness
        EnsCAGskew(n,1)=Features_Ensemble.EnsCAGstats{c,n}(6);
        % Ensemble Co-Activity-Graphy Kurtosis
        EnsCAGkurt(n,1)=Features_Ensemble.EnsCAGstats{c,n}(7);
        % Inter Ensemble Interval Mean ******************
        EnsIEImean(n,1)=Features_Ensemble.IEIstats{c,n}(1);
        % Inter Ensemble Interval Mode 
        EnsIEImode(n,1)=Features_Ensemble.IEIstats{c,n}(2);
        % Inter Ensemble Interval Median
        EnsIEImedian(n,1)=Features_Ensemble.IEIstats{c,n}(3);
        % Inter Ensemble Interval Variance
        EnsIEIvar(n,1)=Features_Ensemble.IEIstats{c,n}(4);
        % Inter Ensemble Interval Skewness
        EnsIEIskew(n,1)=Features_Ensemble.IEIstats{c,n}(5);
        % Inter Ensemble Interval Kurtosis
        EnsIEIkurt(n,1)=Features_Ensemble.IEIstats{c,n}(6);
        % Ensemble Duration Mean ****************************
        EnsEDmean(n,1)=Features_Ensemble.EDstats{c,n}(1);
        % Ensemble Duration Mode
        EnsEDmode(n,1)=Features_Ensemble.EDstats{c,n}(2);
        % Ensemble Duration Median 
        EnsEDmedian(n,1)=Features_Ensemble.EDstats{c,n}(3);
        % Ensemble Duration Variance
        EnsEDvar(n,1)=Features_Ensemble.EDstats{c,n}(4);
        % Ensemble Duration Skewness
        EnsEDskew(n,1)=Features_Ensemble.EDstats{c,n}(5);
        % Ensemble Duration Kurtosis
        EnsEDkurt(n,1)=Features_Ensemble.EDstats{c,n}(6);
    end
    % ********************************************************************
    % Create Table
    TensemblesDetails=table(NeuronsOccupancy,TimeOccupancy,DominanceIndex,...
        EnseRate,...
        EnsCAGauc,EnsCAGmean,EnsCAGmode,EnsCAGmedian,EnsCAGvar,EnsCAGskew,EnsCAGkurt,...
        EnsIEImean,EnsIEImode,EnsIEImedian,EnsIEIvar,EnsIEIskew,EnsIEIkurt,...
        EnsEDmean,EnsEDmode,EnsEDmedian,EnsEDvar,EnsEDskew,EnsEDkurt);
    TensemblesDetails.Properties.VariableNames=HeadersFeaturesEnsembles;
    % Save CSV
    if isdir([FileDirSave,NameDir])
        disp('>>Saving Neural Ensembles Details...')
    else % Create Directory
        disp('>>Directory Created: \Ensemble Features')
        mkdir([FileDirSave,NameDir]);
    end
    writetable(TensemblesDetails,[FileDirSave,NameDir,Experiment,'_',Name,'_byENS.csv'],...
        'Delimiter',',','QuoteStrings',true);
    disp(['Neural Ensemble Details Saved from ',Experiment,'-',Names_Conditions{c}])
end
disp('>> Data Exported @ \Ensemble Features.')
% disp('>>END.')