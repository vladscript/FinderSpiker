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
checkname=1;
while checkname==1
    DefaultPath=[FileDirSave,'Processed Data'];
    if exist(DefaultPath,'dir')==0
        DefaultPath=pwd; % Current Diretory of MATLAB
    end
    [FileName,PathName] = uigetfile('*.mat',[' Pick the Analysis File ',Experiment],...
        'MultiSelect', 'off',DefaultPath);
    dotindex=find(FileName=='.');
    if strcmp(FileName(1:end-4),Experiment)
        checkname=0;
        % SAVE DATA
        disp('>>Updating .mat file...')
        save([PathName,FileName],'Features_Ensemble','Features_Condition',...
            '-append');
        disp([Experiment,'   -> UPDATED (Ensembles Features)'])
    elseif FileName==0
        checkname=0;
        disp('*************DISCARDED************')
    else
        disp('Not the same Experiment!')
        disp('Try again!')
    end
end    
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
                     'RateTrans','RateCycles','SimpleCycles','ClosedCycles','OpenedCycles',...
                     'CAGauc','CAGmean','CAGvar','CAGskew','CAGkurt'};
HeadersFeaturesEnsembles={'Neurons','Time','Domminace','Rate',...
     'ensCAGauc','ensCAGmean','ensCAGvar','ensCAGskew','ensCAGkurt'...
     'IEImean','IEIvar','IEIskew','IEIkurt'...
     'EDmean','EDvar','EDskew','EDkurt'...
     'CAGauc','CAGmean','CAGvar','CAGskew','CAGkurt'};
% ROWS Ensembles are added: Ens1; Ens2; ... Ensj

for c=1:C
    fprintf('>> Creating Table Ensemble Feautures of %s \n',Names_Conditions{c});
    %% General INTEL about Neurla Ensembles
    % Features Columns of the Table ***************************************
    Name=Names_Conditions{c};
    NE=Features_Ensemble.Nenscond(c);                   % N ensembles
    Th=Features_Condition.Threshold;                    % Threshold
    DunnIndx=Features_Condition.Dunn(c);                % Dunns Index
    MaxIntraVec=Features_Condition.MIV(c);              % Max Distance Intra Vectors
    ClassError=Features_Condition.ECV_Cond(c);          % Classification Error
    RateTran=Features_Condition.RateTrans(c);           % Rate Transitions
    RateCycl=Features_Condition.RateCycles(c);          % Rate Cycles
    % [%] Simple Cycles
    Portion_simple=Features_Condition.CyclesType(1,c)/sum(Features_Condition.CyclesType(:,c));
    % [%] Closed Cycles
    Portion_closed=Features_Condition.CyclesType(2,c)/sum(Features_Condition.CyclesType(:,c));
    % [%] Opened Cycles
    Portion_opened=Features_Condition.CyclesType(3,c)/sum(Features_Condition.CyclesType(:,c));
    CAGauc=Features_Condition.CAGstats(c,1);            % CAG AutoCorrelation Coefficient
    CAGmean=Features_Condition.CAGstats(c,2);           % CAG Mean
    CAGvar=Features_Condition.CAGstats(c,3);            % CAG Variance
    CAGskew=Features_Condition.CAGstats(c,4);           % CAG Skewness
    CAGkurt=Features_Condition.CAGstats(c,5);           % CAG Kurtosis
    % ********************************************************************
    % Create Table
    Tensemblesfeatures=table(NE,Th,DunnIndx,MaxIntraVec,ClassError,...
        RateTran,RateCycl,Portion_simple,Portion_closed,Portion_opened,...
        CAGauc,CAGmean,CAGvar,CAGskew,CAGkurt);
    Tensemblesfeatures.Properties.VariableNames=HeadersFeaturesCondition;
    % Save CSV
    if isdir([FileDirSave,NameDir])
        disp('>>Saving...')
        writetable(Tensemblesfeatures,[FileDirSave,NameDir,Experiment,'_',Name,'.csv'],...
            'Delimiter',',','QuoteStrings',true);
        disp(['Saved Ensemble Features: ',Experiment,'-',Names_Conditions{c}])
    else % Create Directory
        disp('>>Directory >Ensemble Features< created')
        disp('>>Saving...')
        mkdir([FileDirSave,NameDir]);
        writetable(Tensemblesfeatures,[FileDirSave,NameDir,Experiment,'_',Name,'.csv'],...
            'Delimiter',',','QuoteStrings',true);
        disp('Ensemble Features Directory Created');
        disp(['Saved Ensemble Features: ',Experiment,'-',Names_Conditions{c}])
    end
    %% Detailed Ensembles INTEL
    % [~,EnseDom]=max(Features_Ensemble.Dominance(c,:));  % Dominant Ensemble
%     for n=1:NE
%         HeadersFeatures{end+1}=['RateEns_',num2str(n)];
%         Tensemblesfeatures(1,end+1)=table(Features_Ensemble.Rate(c,n));
%     end
%     for n=1:NE
%         HeadersFeatures{end+1}=['DomEns_',num2str(n)];
%         Tensemblesfeatures(1,end+1)=table(Features_Ensemble.Dominance(c,n));
%     end
    
end
