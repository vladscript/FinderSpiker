%% Make Data set of Statistics of Network Features
%% 0. Load MAT File
% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end))];
[FileName,PathName] = uigetfile('.mat','Load DATASET of GEPHI Features',...
    CurrentPathOK);
fullFileName = fullfile(PathName, FileName);
load(fullFileName);
% setup
% [NE,NC]=size(GEPHIDATA); % Number of Experiments & Conditions
MoreFeats=true; % enter the loop

%% 1. Chooose Feature from table: clustering, grade, etc, etc
busca=true;
auxn=1;
while busca
    if ~isempty( GEPHIDATA{auxn} )
        NetFeaturesNames=GEPHIDATA{1}.Properties.VariableNames(3:end)';
        busca=false;
    else
        auxn=auxn+1;
    end
end
while MoreFeats    
    [index_var,index_CHECK] = listdlg('PromptString','Select a Network Feature:',...
                        'SelectionMode','single',...
                        'ListString',NetFeaturesNames);
    ActualFeature=NetFeaturesNames(index_var);
    %% 1 Check Empty Tables and Empty-Feature Columns
%     ExsitFeature=zeros(NE,NC);
%     EXPLIST=cell(NE,NC);
%     for e=1:NE
%         for c=1:NC
%             X=GEPHIDATA{e,c};
%             if ~isempty(X)
%                 EXPID=X{1,2};
%                 CurrentFeatures=X.Properties.VariableNames(6:end);
%             else
%                 fprintf('>>Missing Experiment @ %s\n',Names_Conditions{c})
%                 EXPID='no_exp';
%                 CurrentFeatures={'empty'};
%             end
%             EXPLIST{e,c}=EXPID;
%             
%             ColumndIndx=strmatch(ActualFeature,CurrentFeatures);
%             if ~isempty(ColumndIndx)
%                 ExsitFeature(e,c)=1;
%                 fprintf('*')
%             else
%                 fprintf('>>Missing Network parameter: %s @ EXP: %s\n',ActualFeature{1},EXPLIST{e,c});
%             end
%         end
%         fprintf('\n')
%     end
    [ExsitFeature,EXPLIST]=emptycolnettables(GEPHIDATA,Names_Conditions,ActualFeature);
    %% 2. Read data from ALL Tables and ALL Conditions
%     TotalTable=table;
%     ALLpreActive=cell(NE,NC);
%     fprintf('>>Reading data for table: \n')
%     n=1;
%     while sum(ExsitFeature(:))>0
%         % No empty Table
%         if ExsitFeature(n)
%             ActualExp=EXPLIST{n};
%             fprintf('> Experiment: %s in ',ActualExp)
%             [Rows,Cols]=find(strcmp(ActualExp,EXPLIST));
%             preActive=[];
%             for k=1:numel(Cols)
%                 ActualCondition=Names_Conditions{Cols(k)};
%                 ActualActive=[];
%                 fprintf('%s,',ActualCondition)
%                 X=GEPHIDATA{Rows(k),Cols(k)};     % Table
%                 x=X{:,{ActualFeature{1}}};  % Cell Values
%                 y=X{:,'label'};     % labels bleonging to ensemble
%                 % Make it Vector: (yeahp, with a loop)
%                 xnum=zeros(numel(x),1);
%                 for nx=1:numel(x)
%                     xnum(nx)=str2double(x{nx});
%                     if y{nx}~='0'
%                         ActualActive=[ActualActive;nx];
%                     end
%                 end
%                 % Only Previus and Actual Nodes #######################
%                 OKindex=union(preActive,ActualActive);      % join
%                 preActive=OKindex;                  % update
%                 RowTable=table({ActualCondition},{ActualExp},...
%                     mean(xnum(OKindex)),mode(xnum(OKindex)),...
%                     median(xnum(OKindex)),var(xnum(OKindex)),...
%                     skewness(xnum(OKindex)),kurtosis(xnum(OKindex)));
%                 % Statistics: mean, variance, mode, median, etc
%                 TotalTable=[TotalTable;RowTable];
%                 ExsitFeature(Rows(k),Cols(k))=0;
%             end
%             fprintf('\n')
%         end
%         n=n+1;
%     end
%     fprintf('>>Data Loaded.\n')
    [TotalTable,DATAnet]=TableNetwork(ExsitFeature,GEPHIDATA,EXPLIST,Names_Conditions,ActualFeature);
    
%% 5. Display Data Statistics
            % Set Name to Variables
    TotalTable.Properties.VariableNames={'Condition','EXPID',...
        ['mean_',ActualFeature{1}],['mode_',ActualFeature{1}],...
        ['median_',ActualFeature{1}],['var_',ActualFeature{1}],...
        ['skew_',ActualFeature{1}],['kurt_',ActualFeature{1}]};
% Plot @ >>Feature_Explorer
% %     % Plot Figure
% %     FigNetStats=figure;
% %     FigNetStats.Name=['Network ',ActualFeature{1},' Stats'];
% %     FigNetStats.Position=[55,207,891,420];
% %     titlesSubs={'mean','mode','median','var','skew','kurt'};
% %     % COLORS MAP
% %     SetColorMap; % SAME AS AIMs colors
% %     ColorsMap=cbrewer(KindMap,ColorMapName,Ncolors);
% %     % Condition List:
% %     T1=TotalTable{:,'Condition'};
% %     IndexesCond=zeros(NE,NC);
% %     DataMatrix=zeros(NE,NC);
% %     hlabel=[];
% %     for nstat=1:6;
% %         step=-0.05;
% %         for c=1:NC
% %             IndexesCond(:,c)=strmatch(Names_Conditions{c},T1);
% %             DataMatrix(:,c)=TotalTable{IndexesCond(:,c),2+nstat};
% %             % PLOT RAIN CLOUDS
% %             subplot(3,2,nstat)
% %             step=step+.2;
% %             hplot{nstat,c}=raincloud_plot(DataMatrix(:,c),'color',ColorsMap(c,:),'box_on',1,'alphaval',0.5,...
% %              'box_dodge', 1, 'box_dodge_amount',step , 'dot_dodge_amount', step, 'box_col_match',1,...
% %              'band_width',0.2);
% %             axis tight; grid on;
% %             title(titlesSubs{nstat});
% %             if nstat==1
% %                 % save axis object to legend
% %                 hlabel=[hlabel,hplot{nstat,c}{1}];
% %             end
% %         end
% %     end
% %     legend(hlabel,Names_Conditions,'Location','northwestoutside');
    %% 6. Make Table for ML o Statistical Anlysis
    okbutton = questdlg('Make CSV Table?');
    waitfor(okbutton); 
    if strcmp('Yes',okbutton)
        % Netwrok Directory
        CurrentPath=pwd;
        Slshes=find(CurrentPath=='\');
        CurrentPathOK=[Dirpwd(1:slashesindx(end)),'NetWorks-CSV'];
        % Set Save Name
        timesave=clock;
        TS=num2str(timesave(1:5));
        TS=TS(TS~=' ');
        SaveFile=['\Network_',ActualFeature{1},'_Dataset_',TS,'.csv'];
        % Select Destiny
        PathSave=uigetdir(CurrentPathOK);
        disp('>>Making CSV table...')
        writetable(TotalTable,[PathSave,SaveFile],...
                        'Delimiter',',','QuoteStrings',true);
        fprintf('>> Dataset saved @: %s\n',[PathSave,SaveFile])
    else
        fprintf('>>Unsaved dataset.\n')
    end

    answer = questdlg('Inspect More Features?', ...
    'NETWROK FEATURES', ...
    'YES','NO','NO');

    switch answer
        case 'YES'
            disp([answer ' more Network features coming right up.'])
            MoreFeats = true;
        case 'NO'
            disp([answer ' more Network features.'])
            MoreFeats = false;
    end

    % Ask if Save CSV, clear and Goodbye
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')
%% END OF THE WORLD