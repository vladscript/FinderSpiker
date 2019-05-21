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
[NE,NC]=size(GEPHIDATA); % Number of Experiments & Conditions
MoreFeats=true; % enter the loop

%% 1. Chooose Feature from table: clustering, grade, etc, etc
NetFeaturesNames=GEPHIDATA{1}.Properties.VariableNames(3:end)';
while MoreFeats    
    [index_var,index_CHECK] = listdlg('PromptString','Select a Network Feature:',...
                        'SelectionMode','single',...
                        'ListString',NetFeaturesNames);
    ActualFeature=NetFeaturesNames(index_var);
    %% 1.1 Check if the features is in ALL tables
    ExsitFeature=zeros(NE,NC);
    EXPLIST=cell(NE,NC);
    for e=1:NE
        for c=1:NC
            X=GEPHIDATA{e,c};
            if ~isempty(X)
                EXPID=X{1,2};
                CurrentFeatures=X.Properties.VariableNames(6:end);
            else
                fprintf('>>Missing Experiment @ %s\n',Names_Conditions{c})
                EXPID='no_exp';
                CurrentFeatures={'empty'};
            end
            EXPLIST{e,c}=EXPID;
            
            ColumndIndx=strmatch(ActualFeature,CurrentFeatures);
            if ~isempty(ColumndIndx)
                ExsitFeature(e,c)=1;
                fprintf('*')
            else
                fprintf('\n>>Missing Feature: %s @ EXP: %s\n',ActualFeature{1},EXPLIST{e,c});
            end
        end
        fprintf('\\')
    end
    if sum(ExsitFeature(:))==NE*NC
        fprintf('\n>>Data Correctly Read for <%s> Network Feature\n',ActualFeature{1})
        goon=true;
    else
        disp('>>Missing Data, please re-check your tables')
        goon=false;
    end
    if goon
        %% 2. Check Experiment ID's are OK (non repeated and matched)
        AreExps=zeros(NE,1);
        IndexesExps=zeros(NE,NC);
        IndexesExps(:,1)=1:NE;
        counterok=0;
        for e=1:NE
            FirstList=EXPLIST{e,1};
            for c=2:NC
                NextList=EXPLIST(:,c);
                IndxExp=strmatch(FirstList,NextList);
                if ~isempty(IndxExp)
                    IndexesExps(e,c)=IndxExp;
                    counterok=counterok+1;
                else
                    disp('>>Error at list of Experiments')
                end
            end
            fprintf('*\\')
        end
        fprintf('\n>>List Of Experiments: Alright!')
        if counterok==NE*(NC-1)
            gomoreon=true;
        else
            gomoreon=false;
        end
        if gomoreon
            %% 3. Read vector from ALL Tables and ALL Conditions
            TotalTable=table;
            fprintf('>>Loading: \n')
            for e=1:NE
                ActualExp=EXPLIST{e,1};
                for c=1:NC
                    ActualCondition=Names_Conditions{c};
                    RowIndx=IndexesExps(e,c);   % Index of the Experiment
                    X=GEPHIDATA{RowIndx,c};     % Table
                    x=X{:,{ActualFeature{1}}};  % Cell
                    % Make it Vector: (yeahp, with a loop)
                    xnum=zeros(numel(x),1);
                    for nx=1:numel(x)
                        xnum(nx)=str2double(x{nx});
                    end
                    RowTable=table({ActualCondition},{ActualExp},...
                        mean(xnum),mode(xnum),median(xnum),var(xnum),skewness(xnum),kurtosis(xnum));
                    %% 4. Make Statistics: mean, variance, mode, median, etc
                    TotalTable=[TotalTable;RowTable];
                    fprintf('*')
                end
                fprintf('\\')
            end
            disp('>>Loaded.')
            %% 5. Display Data Statistics
            % Set Name to Variables
            TotalTable.Properties.VariableNames={'Condition','EXPID',...
                ['mean_',ActualFeature{1}],['mode_',ActualFeature{1}],...
                ['median_',ActualFeature{1}],['var_',ActualFeature{1}],...
                ['skew_',ActualFeature{1}],['kurt_',ActualFeature{1}]};
            % Plot Figure
            FigNetStats=figure;
            FigNetStats.Name=['Network ',ActualFeature{1},' Stats'];
            FigNetStats.Position=[55,207,891,420];
            titlesSubs={'mean','mode','median','var','skew','kurt'};
            % COLORS MAP
            SetColorMap; % SAME AS AIMs colors
            ColorsMap=cbrewer(KindMap,ColorMapName,Ncolors);
            % Condition List:
            T1=TotalTable{:,'Condition'};
            IndexesCond=zeros(NE,NC);
            DataMatrix=zeros(NE,NC);
            hlabel=[];
            for nstat=1:6;
                step=-0.05;
                for c=1:NC
                    IndexesCond(:,c)=strmatch(Names_Conditions{c},T1);
                    DataMatrix(:,c)=TotalTable{IndexesCond(:,c),2+nstat};
                    % PLOT RAIN CLOUDS
                    subplot(3,2,nstat)
                    step=step+.2;
                    hplot{nstat,c}=raincloud_plot(DataMatrix(:,c),'color',ColorsMap(c,:),'box_on',1,'alphaval',0.5,...
                     'box_dodge', 1, 'box_dodge_amount',step , 'dot_dodge_amount', step, 'box_col_match',1,...
                     'band_width',0.2);
                    axis tight; grid on;
                    title(titlesSubs{nstat});
                    if nstat==1
                        % save axis object to legend
                        hlabel=[hlabel,hplot{nstat,c}{1}];
                    end
                end
            end
            legend(hlabel,Names_Conditions,'Location','northwestoutside');
            %% 6. Make Table for ML o Statistical Anlysis
            okbutton = questdlg('Make CSV Table?');
            waitfor(okbutton); 
            if strcmp('Yes',okbutton)
                % Netwrok Directory
                CurrentPath=pwd;
                Slshes=find(CurrentPath=='\');
                CurrentPathOK=[CurrentPath(1:Slshes(end)),'NetWorks-CSV'];
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
    end
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')
%% END OF THE WORLD