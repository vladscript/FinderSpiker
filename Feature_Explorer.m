%% Feature Explorer ********************************
% Script to plot boxplots, pdfs and explore effects
% 
%% Load Tabel **************************************
Import_FinderSpiker;
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end))]; % Finder Spiker Main Folder
% Load File 
[FileName,PathName,MoreFiles] = uigetfile({'*.csv'},' Dataset file of ALL Features',...
    'MultiSelect', 'off',CurrentPathOK);

Xraw=readtable([PathName,FileName]);                % Table
%% Setup :
Y=categorical(table2array( Xraw(:,1)) );            % Labels
EXPIDs=table2array( Xraw(:,2));                     % Cell of Strings
X=table2array( Xraw(:,3:end) );                     % Dataset
FeatureNames=Xraw.Properties.VariableNames(3:end);  % Feature Names
Labels=unique(Y);                                   % Conditions Names
Nconditions=numel(Labels);                          % N conditons
% Set Names Coniditons:
for c=1:Nconditions
    NamesCond{c}=char(Labels(c));                   % Cell arrray of Strings
end
Nfeat=numel(FeatureNames);                          % N features
%% Color Selector
% Colors: #################################################################
[CM,IndxColor]=Color_Selector(NamesCond);

%% p-values: descriptive statistical tests

% Set Statistical Tes
if Nconditions<3
    % Two Paired Groups: 
    % Rank, Score, or Measurement (from Non- Gaussian Population)
    Test2Do='Paired';
else
    % More than 2 Conditions
    % Rank, Score, or Measurement (from Non- Gaussian Population)
    Test2Do='Repeated';
end

%% FEATURE EXPLORER #######################################################
% Save Results to Retrieve Them
DataCell=cell(Nconditions,1);
for f=1:Nfeat % Feature Loop **********************************************
    figure;
    labelsplot=[];
    DodgeInit=0.3;
    Dodge=DodgeInit;
    % Set Text, Retrieve & Plot Data
    data=X(:,f);
    fprintf('>> %s @ Conditions:\n',char(FeatureNames(f)));
    for c=1:Nconditions
        if c<Nconditions
            fprintf('  %s vs',char(Labels(c)))
        else
            fprintf('  %s',char(Labels(c)))
        end
        DataCell{c}=data(Y==Labels(c));
        hplot{c}=raincloud_plot(DataCell{c},'color',CM(IndxColor(c),:),'box_on',1,'alphaval',1,'box_dodge', 1, 'box_dodge_amount',Dodge, 'dot_dodge_amount', Dodge, 'box_col_match',0);
        labelsplot=[labelsplot,hplot{c}{1}];
        Dodge=DodgeInit*(c+1);
        LabelName{c}=char(Labels(c));
    end
    axis tight; grid on;
    legend(labelsplot,LabelName);
    title(char(FeatureNames(f)))
    fprintf('\n')
    
    % Statistics Data
    
%     switch Test2Do
%         case 'Paired'
%             % 2 Conditions
%             A=DataCell{1};
%             B=DataCell{2};
%             [p,h] = ranksum(A,B);
%         case 'Repeated'
%             % Make CODE!
%             datastat=reshape(cell2mat(DataCell),[size(DataCell{1},1),size(DataCell,1)]);
%             [p,tbl,stats] = friedman(datastat,size(DataCell,1));
%             c=multcompare(stats);
%     end
%     fprintf('>>P-value: p=%f\n',p);
    pause;
end
%% For paired Experiments: Delta Features
% find and compare paired Features
% Loop for unique Expeiment IDs: if more than 1-> paired experiment

Labels=unique(Y);
PairedExps=cell(numel(Labels));
for lrow=1:numel(Labels)
    for lcol=lrow+1:numel(Labels)
        % List of Experiments for Condition A
        ActualLabel=Labels(lrow);
        IndxConditionA=find(Y==ActualLabel);
        EXPlistA=unique(EXPIDs(IndxConditionA));
        % List of Experiments for Condition B
        NextLabel=Labels(lcol);
        IndxConditionB=find(Y==NextLabel);
        EXPlistB=unique(EXPIDs(IndxConditionB));
        % Intersection of Experiments for Condition A & B
        PairedAB=intersect(EXPlistA,EXPlistB);
        if ~isempty(PairedAB)
            fprintf('>>Paired Experiments for %s & %s:\n',char(ActualLabel),...
                char(NextLabel))
            disp(PairedAB)
        else
            fprintf('>>No Paired Experiments for %s & %s:\n',char(ActualLabel),...
                char(NextLabel))
        end
        PairedExps{lrow,lcol}=PairedAB;
    end
end

% % %% Calculate Deltas ******************************************************
% % DeltaExps=cell(numel(Labels));
% % ReferenceCondition={};
% % for lrow=1:numel(Labels)
% %     for lcol=lrow+1:numel(Labels)
% %         if ~isempty(PairedExps{lrow,lcol})
% %             EXPlist=PairedExps{lrow,lcol};
% %             LabelsDelta={char(Labels(lcol));char(Labels(lrow))};
% %             % Select Conditions in Order To Calculate Deltas
% %             for c=1:2
% %                 [index_var(c),~] = listdlg('PromptString',...
% %                     ['Set Condition in Order: ',num2str(c)],...
% %                     'SelectionMode','single',...
% %                     'ListString',LabelsDelta);
% %             end
% %             DeltaCondition=LabelsDelta(index_var);
% %             % Save Reference Condition
% %             ReferenceCondition=[ReferenceCondition;DeltaCondition(1)]
% %             % For every Experiment:
% %             DeltaAllFeature=[];
% %             for e=1:numel(EXPlist)
% %                 IndxTable=find(ismember(EXPIDs,EXPlist{e}));
% %                 % Cond_A - Cond_B: COND_B/COND_A
% %                 RowA=IndxTable(Y(IndxTable)==DeltaCondition(1));
% %                 RowB=IndxTable(Y(IndxTable)==DeltaCondition(2));
% %                 % Relative Changes
% %                 DeltaFeature(X(RowB,:)~=0)=X(RowB,X(RowB,:)~=0)./X(RowA,X(RowB,:)~=0);
% %                 DeltaFeature(X(RowB,:)==0)=X(RowA,X(RowB,:)==0)-X(RowB,X(RowB,:)==0);
% %                 DeltaAllFeature=[DeltaAllFeature;100*DeltaFeature];
% %             end
% %             DeltaExps{lrow,lcol}=DeltaAllFeature;
% %         else
% %             disp('No Paired Conditions')
% %         end
% %     end
% % end
% % ReferenceCondition=unique(ReferenceCondition); % Condition as References
% % %% Show Results
% % RefIndex=[];
% % for c=1:numel(ReferenceCondition)
% %     RefIndex=[RefIndex;find(Labels==ReferenceCondition{c})];
% % end
% % 
% % for c=1:numel(ReferenceCondition)
% %     TitleFig=Labels(RefIndex(c));
% %     DeltaNum=[]; ConditionLabel={}; % To Make Boxplots
% %     for r=1:numel(Labels)
% %         if ~isempty(DeltaExps{r,RefIndex(c)})
% %            % Gather All The Deltas 
% %            VersusCondition=Labels(r);
% %            Nexps=size(DeltaExps{r,RefIndex(c)},1);
% %            for e=1:Nexps
% %                 ConditionLabel=[ConditionLabel;['+ ',char(VersusCondition)]];
% %            end
% %            DeltaNum=[DeltaNum;DeltaExps{r,RefIndex(c)}];
% %         end
% %     end
% %     figure;
% %     for IndexFeat=1:size(DeltaNum,2)
% %         boxplot(DeltaNum(:,IndexFeat),ConditionLabel);
% %         ylab=ylabel(['%\Delta',FeatureNames{IndexFeat}]);
% %         ylab.Interpreter='tex';
% %         title(['Reference: ',char(TitleFig)])
% %         pause;
% %     end
% % end
