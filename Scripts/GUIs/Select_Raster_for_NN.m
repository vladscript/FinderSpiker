% Function to Select frames per Condition Manually
% Extract Descriptive RasterFeatures and seve the, @
% ...\FinderSpiker\Raster Features
% Run After Finder_Spiker_Calcium.m and Detected & Undetected Scripts
% Load Saved Data in the EXPERIMENT_DATE.mat Files
% It runs until selection is OK (according to the user)
% Input
%   fs:                 Sampling Frequency
%   Raster_Condition:   Separated Raster in Cells   NOT condition SORTED
%   XY:                 Original Coordinates        NOT condition SORTED
%   Names_Conditions:   Names of Conditions
%   Experiment:         Name of the Experiment
%   ESTSIGNAL:          Denoised Signals (transients counting)
%   checkname           if it's specified->OMITS SAVING
% Ouput
%   RASTER_Selected_Clean   Matrix (ONLY ACTIVE cells)
%   XY_selected:            Selected ACtive Cells
%   R_Condition:            Cell of Raster for each COndition
%   Onsets:                 Selected Starting Point
%   New_indexes:            Sleected & Sorted Indexes of the Cells
function [RASTER_Selected_Clean,XY_selected,R_Condition,Onsets,New_Index]= Select_Raster_for_NN(fs,Raster_Condition,XY,Names_Conditions,Experiment,ESTSIGNALS,varargin)
%% Setup
dontupdmat=false;
if numel(varargin)==1
    if ismatrix(varargin{1})
        checkname=1;
        Onsets=varargin{1};
        RasterDurations=get_raster_durations(Onsets,Raster_Condition,fs);
        dontupdmat=true;
    else
        checkname=0;
    end
        
else
    checkname=1;
end
FileDirSave=pwd;
Experiment=Experiment(Experiment~='\'); % PATCH
slashes=find(FileDirSave=='\');
FileDirSave=FileDirSave(1:slashes(end));
okbutton='No';
while ~strcmp('Yes',okbutton)
     %% Show Raster per condition
    [~,NC]=size(Raster_Condition);
    Raster_Selecter={};
    Onsets=cell(NC,1);
    XDEN=cell(NC,1); % Denoised Signals per Condition
    Ncells=size(ESTSIGNALS{1},1); % NUmber of Cells 
    for c=1:NC
        R=Raster_Condition{c};
        XdenoisedAll=[]; % Denoised Signals Buffer
        [~,FramesSize]=size(R);
        Plot_Raster_Ensembles(R,fs);
        figure_raster=gcf;
        figure_raster.Name=['Raster ',Names_Conditions{c}];
        % Read Selection
        prompt = {'Start @ min:','Finish @ min:'};
        dlg_title = 'Select Raster ';
        num_lines = 1;
        defaultans = {'0',num2str(ceil(FramesSize/fs/60))};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Minute_Start=str2double(answer{1});
        Minute_End=str2double(answer{2});
        Frame_A=round(Minute_Start*60*fs+1);
        Frame_B=round(Minute_End*60*fs);
        Onsets{c,1}=Frame_A;
        if Frame_B>FramesSize
            Frame_B=FramesSize;
        end
        % Selection:
        R=R(:,Frame_A:Frame_B);
        Raster_Selecter{c}=R; 
        % Denoised Signal
        for j=1:size(ESTSIGNALS,1)
            XdenoisedAll=[XdenoisedAll,ESTSIGNALS{j,c}];
        end
        if size(Raster_Condition{c},2)~=size(XdenoisedAll,2)
            Frame_A=round(RasterDurations(c,1)*60*fs+1);
            Frame_B=round(RasterDurations(c,2)*60*fs);
        end
        XDEN{c}=XdenoisedAll(:,Frame_A:Frame_B);
        % Show Selected Raster
        close(figure_raster);
        NVaux{c,1}='1'; % Number of "Videos" or Records
    end
    %% Ignore NonActive Coordinates in the Selection and Sort per Condition
    % Sort Neurons by Activation in Each Condition:
    [New_Index,Raster_Condition_Sel,RASTER_WHOLE]=SortNeuronsCondition(Raster_Selecter,Ncells);
    % Plot Whole Sorted Raster:
    
    XY_selected=XY(New_Index,:);
    % Check repeated coordinates
    IndxRepeated=repxychekcker(XY_selected);
    if ~isempty(IndxRepeated)
        fprintf('\n>Discarding repeated coordinates: ')
        New_Index=New_Index(setdiff(1:size(XY_selected,1),IndxRepeated(2:2:end)));
    end
    fprintf('ready\n')
    XY_selected=XY(New_Index,:);
    RASTER_Selected_Clean=RASTER_WHOLE(New_Index,:);
    % Clean Raster and Coordinates
    % INDEX of SORTED Active NEURONS during the WHOLE experiment
    ActiveNeurons=find(sum(RASTER_Selected_Clean,2)>0);                 
    New_Index(ActiveNeurons)
    RASTER_Selected_Clean=RASTER_Selected_Clean(ActiveNeurons,:);
    XY_selected=XY_selected(ActiveNeurons,:);
    Plot_Raster_Ensembles(RASTER_Selected_Clean,fs);
    set(gcf,'Name',['ID: ',Experiment,' selected '],'NumberTitle','off')
    Label_Condition_Raster(Names_Conditions,Raster_Condition_Sel,fs);
    % To Save active and Sorted:
    R_Condition={};

    %% Feature Table Column Names
    HeadersFeatures={'RateNeurons','ActivityTimeFraction','ActiveRatioCAG','EffectiveActivity',...
        'ISImean','ISImode','ISImedian','ISIvar','ISIskew','ISIkurt',...
        'Lengthmean','Lengthmode','Lengthmedian','Lengthvar','Lengthskew','Lengthkurt',...
        'CAGmean','CAGmode','CAGmedian','CAGvar','CAGskew','CAGkurt',...
        'RoAmean','RoAmode','RoAmedian','RoAvar','RoAskew','RoAkurt',...
        'RoTmean','RoTmode','RoTmedian','RoTvar','RoTskew','RoTkurt'};
    %% FEATURES EXTRACTION and plot
    preActive=[];
    for c=1:NC
        % Read Original Raster (UNSORTED!!!)
        R=Raster_Condition_Sel{c};
        % Raster (SORTED!!!)
        R=R(New_Index,:);
        Xden=XDEN{c}(New_Index,:);
        % Raster (JUST ACTIVE of the WHOLE EXPERIMENT!!!)
        R_Condition{c}=R(ActiveNeurons,:);
        Xden=Xden(ActiveNeurons,:);
        % XY_Condition{c}=XY_selected;
        % Only Actual And Previous
        ActualActive=find(sum(R_Condition{c},2)./size(R_Condition{c},2));     % read
        OKindex=union(preActive,ActualActive);      % join
        preActive=OKindex;                          % update
        % Activity Indexes:*******************************************
        % Feature of ALL Rows of the SELECTED RASTER (BETTER OFF!)
        [Descriptive,AI_1(c,1),AI_2(c,1)]=get_general_features(R_Condition{c}(OKindex,:));
        AN=Descriptive.AN;       % Active Neurons
        DurAc=Descriptive.DurAc; % Number of Active Frames
        CAG=Descriptive.CAG;     % Coactivity
        RoA=Descriptive.RoA;     % Active Frames / Total Frames
        RAN=Descriptive.RAN;     % Rate of Active Neurons
        
        % delete:__________________________________________________________
        % [~,~,~,~,StatsFeatures]=get_iti_pdf(R_Condition{c},fs);
        % _________________________________________________________________
        % Temporal Measures According to Denoised Signal (wavelet analysis)
        [NoT,ITI,LT]=get_RoT(R_Condition{c}(OKindex,:),Xden);
        ITI=ITI/fs; LT=LT/fs;                                   % [SECONDS]
        RasterDuration=size(R_Condition{c},2)/fs/60;            % [MINUTES]
        RoT=NoT/RasterDuration; % Ca++Tranisetns per minute
        % -----------------------------------------------------------------
        % CoActivityGram Statistics & pdf: IMPORTANT!!!!
%         if or(max(CAG)==0,min(CAG)==max(CAG))
%             CAGbin=linspace(min(CAG),max(CAG),100);
%             CAGp=zeros(size(CAGbin));
%         else
%             [~,~]=ksdensity(CAG,linspace(min(CAG),max(CAG),100));
%         end
        % Feature Table *******************************************
        % 'RateofNeurons','ActivityDuration','MeanActivity','EffectiveActivity',...
        % 'ISImean','ISImode','ISImedian','ISIvar','ISIskew','ISIkurt',...
        % 'Lengthmean','Lengthmode','Lengthmedian','Lengthvar','Lengthskew','Lengthkurt',...
        % 'CAGmean','CAGmode','CAGmedian','CAGvar','CAGskew','CAGkurt',...
        % 'RoAmean','RoAmode','RoAmedian','RoAvar','RoAskew','RoAkurt',...
        % 'RoTmean','RoTmode','RoTmedian','RoTvar','RoTskew','RoTkurt']
        Trasterfeatures=table(RAN,DurAc,AI_1(c,1),AI_2(c,1),...
            mean(ITI),mode(ITI),median(ITI),var(ITI),skewness(ITI),kurtosis(ITI),...
            mean(LT), mode(LT), median(LT), var(LT), skewness(LT), kurtosis(LT),...
            mean(CAG),mode(CAG),median(CAG),var(CAG),skewness(CAG),kurtosis(CAG),...
            mean(RoA),mode(RoA),median(RoA),var(RoA),skewness(RoA),kurtosis(RoA),...
            mean(RoT),mode(RoT),median(RoT),var(RoT),skewness(RoT),kurtosis(RoT),...
            'VariableNames',HeadersFeatures);
        disp('>>Raster Features Table: Done.');
        % Saving CSV - - - - - - - - - - - - - - - - - - - - - - - - -
        if checkname
            NameDir='Raster Features\';
            Experiment=Experiment(Experiment~='\');
            if isdir([FileDirSave,'\Raster Features'])
                disp(['Saved Raster Features: ',Experiment,'-',Names_Conditions{c}])
            else % Create Directory
                disp('Directory >Raster Features< created')
                mkdir([FileDirSave,NameDir]);
                
                disp('Resume Tables Directory Created');
                disp(['Saved Raster Features: ',Experiment,'-',Names_Conditions{c}])
            end
            writetable(Trasterfeatures,[FileDirSave,NameDir,Experiment,'_',Names_Conditions{c},'Raster_Features.csv'],...
                    'Delimiter',',','QuoteStrings',true);
        end
    end
    okbutton = questdlg('Selection Alright?');
    %% SAVE OUTPUT DATASET (.m file)
    %checkname=1; % USE AS INPUT
    while ~dontupdmat && strcmp('Yes',okbutton)
        DefaultPath=pwd; % Current Diretory of FinderSpiker
        slashes=find(DefaultPath=='\');
        DefaultPath=[DefaultPath(1:slashes(end)),'\Processed Data'];        
        [FileName,PathName] = uigetfile('*.mat',[' Pick the Analysis File to SAVE CHANGES of: ',Experiment],...
            'MultiSelect', 'off',DefaultPath);
        % dotindex=find(FileName=='.');
        if strcmp(FileName(1:end-4),Experiment)
            checkname=0;
            % SAVE DATA
            save([PathName,FileName],'RASTER_Selected_Clean','XY_selected',...
                'R_Condition','New_Index','Onsets','-append');
            disp([Experiment,'   -> UPDATED (Selected Data)'])
            dontupdmat=true;
        elseif FileName==0
            checkname=0; dontupdmat=true;
            disp('*************DISCARDED************')
        else
            disp('Not the same Experiment!')
            disp('Try again!')
            % okbutton='No';
        end
    end    
end
%% END OF THE WORLD