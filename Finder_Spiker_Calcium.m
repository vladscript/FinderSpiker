% Main script to Load & Process Calcium Fluorescence Signals
% It detects Calcium Transients from VIDEOS of a single slice (EXPERIMENT)
% Sugggested Directorie's Strcuture:
% NAME_EXPERIMENT/ {VIDEOS & Coordinates from fSIENN):
% Input (same directory): 
%   ALL (XY).csv        Coordinates from 4th row x,y,r
%   Condition_A_00.avi
%   Condition_A_01.avi
%   ...
%   Condition_B_00.avi
%   ... 
%   Condition_Z_##.avi
% Output
%   ExperimentID.mat @ cd ..\Processed Data: Useful Workspace Variables
%   ExperimentID.csv @ cd ..\Features Tables:Processing Signal Features
%   ExperimentID.csv @ cd ..\Resume Tables:  Experiment Resume Features
%   Using GitHub                    19/01/2018
%   Important Update:               09/07/2018
%% Global Setup ***********************************************************
clc
clear;
close all;
%% Global Variables:
global SIGNALS;
global DETSIGNALS;
global preDRIVE;
global preLAMBDAS;
global RASTER
global Responses
global Names_Conditions;
global fs;
global SIGNALSclean;
global SNRlambda;
global Experiment;
global RasterAlgorithm;
% global isSIGNALS;
% global notSIGNALS;
%% ADDING ALLSCRIPTS
Update_Directory

%% Set Default Directory of Experiments
DefaultPath='C:\Users\Vladimir\Documents\Doctorado\Experimentos\'; % Load from DEFAULT
if exist(DefaultPath,'dir')==0
    DefaultPath=pwd;
end
%% Read Sampling Frequency
fs=NaN; % To read fs to get into the following loop:
while isnan(fs) % Interface error user reading fs
    fs = inputdlg('Sampling Frequency [Hz] : ',...
                 'VIDEOS', [1 50]);
    fs = str2double(fs{:});
end
% Read Fluorophore DYe
dyename = inputdlg('Fluorophore : ',...
             'DYE', [1 50]);
% fs = str2double(fs{:});

%% Read Names, Path and Coordinates***********************************
[Names_Conditions,NumberofVideos,FN,PathName,XY,r]=Read_Videos(DefaultPath);
for v=1:length(NumberofVideos)
    NVal(v)=round(str2double(NumberofVideos(v)));
end
NV=max(NVal);   % Max N of Videos
[~,NC]=size(FN);                                % N Conditions
%% Initalization Data Output
SIGNALS=cell(NV,NC);
DETSIGNALS=cell(NV,NC);
ESTSIGNALS=cell(NV,NC);
SNRwavelet=cell(NV,NC);
Responses=cell(NV,NC);
preDRIVE=cell(NV,NC);
preLAMBDAS=cell(NV,NC);
TAUSall=cell(NV,NC);
% SIGNALSclean=cell(NV,NC);
RASTER=cell(NV,NC);
isSIGNALS=cell(NV,NC);
notSIGNALS=cell(NV,NC);
LAMBDASSpro=cell(NV,NC);
SNRs=cell(NV,NC);
DRIVERs=cell(NV,NC);
SIGNALSclean=cell(size(SIGNALS)); % Detected clean SIGNALS
SNRlambda=cell(size(preLAMBDAS));    % Sparse Empirical SNRs
RasterAlgorithm='Driver'; % 'OOPSI', 'Derivative'
%% Load Data *********************************************************
% For each CONDITION and VIDEO
for i=1:NC
    for j=1:str2double(NumberofVideos{i})
        FileName=FN{j,i};
        [mov]=Video_Load(FileName,PathName);    % Load Video
        [FS]=Fluorescence_Load(mov,XY,r);       % Load Fluorescence
        SIGNALS{j,i}=FS;
    end
    disp('***')
end
[H,W]=size(mov(1).cdata);   % Height & Width
clear mov;                  % Clear Video Structure
%% Save(1) RAW Data * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% Direcotry to Save: Up from  this one (pwd)
FileDirSave=pwd;
slashes=find(FileDirSave=='\');
FileDirSave=FileDirSave(1:slashes(end));
%% Save data
% Get the Experiment ID:
slashes=find(PathName=='\');
Experiment=PathName(slashes(end-1):slashes(end)-1); % Experiment ID
if isdir([FileDirSave,'\Processed Data'])
    save([FileDirSave,'\Processed Data',Experiment,'.mat'],'Experiment','SIGNALS',...
    'Names_Conditions','NumberofVideos','XY','fs','r');
    disp('SAVED RAW DATA')
else % Create Directory
    disp('Directory >Processed Data< created')
    mkdir([FileDirSave,'\Processed Data']);
    save([FileDirSave,'\Processed Data',Experiment,'.mat'],'Experiment','SIGNALS',...
    'Names_Conditions','NumberofVideos','XY','fs','r');
    disp('SAVED RAW DATA')
end

%% SETUP PROCESSING PARAMTERS ******************************
% Setup for Auto Regressive Process Estimation
% L=30;              % seconds  of the Fluorophore Response
% p=3;               % AR(p) initial AR order        
% taus_0= [.75,2,1]; % starting values for taus
Load_Default_Values_SP;

%% DETRENDING AND SPARSE DECONVOLVE [AUTOMATIC] DETECTION
[NV,NC]=size(SIGNALS);
ColumnNames={'Fluo_Dye','f_s','DetectedCells','Frames',...
        'minSNR','minSkewness','TimeProcessing'};
disp('**************** +++[Processing]+++ *************************')
for i=1:NC
% for i=2:NC
    auxj=0;
    for j=1:str2double(NumberofVideos{i})
        % Initialize and Read Data
        X=SIGNALS{j,i};
        [Ns,Frames] = size(X);
        FR=zeros(Ns,round(L*fs));
        TAUS=zeros(Ns,3);
        DRIVER=zeros(Ns,Frames);
        LAMBDASS=zeros(Ns,1);
        Dfix=zeros(Ns,Frames);
        % XDfix=zeros(Ns,Frames);
        Xestfix=zeros(Ns,Frames);
        XDupdate=zeros(Ns,Frames);
        % Display Info **********************
        disp(['               [>> Condition: ',num2str(i),'/',num2str(NC),']']);
        disp(['               [>> Video    : ',num2str(j),'/',NumberofVideos{i},']']);
        tic
        % Detrending ******************************************************
        % XD=only_detrending(X);
        % Denoising & Feature Extration ***********************************
        [Xest,SNRbyWT,SkewSignal,~,SkewNoise,XDupdate]=denoise_wavelet(X);
        % Find Basal Slow Changing Signals
        % BasalSwitch=basaldetector(X,XDupdate,Xest);
        %% DRIVER ANALYSIS
        IndexesFix=1;       % TO enter to the WhileLoop
        FixedSignals=[];    % 
        FixedNOSignals=[];  % 
        aux=0;              % Loop COunter 
        while and(~isempty(IndexesFix),aux<2)
            aux=aux+1;
            if ~isempty(FixedSignals);
                SkewSignal(FixedSignals)=Features(:,1);
                SkewNoise(FixedSignals)=Features(:,2);
                SNRbyWT(FixedSignals)=Features(:,3);
            end
            if ~isempty(FixedNOSignals)
                SkewSignal(FixedNOSignals)=FeaturesR(:,1);
                SkewNoise(FixedNOSignals)=FeaturesR(:,2);
                SNRbyWT(FixedNOSignals)=FeaturesR(:,3);
            end
            %%% Decision Features
            % SNR >0 *******************************************
            SNRindx=find(SNRbyWT>0);
            % Skew PDFs ****************************************
            % [skew'noise PDF is VERY like RANDN(Ns;Frames)]
            Th_Skew=max(SkewNoise);
            indxSKEW=find(SkewSignal>Th_Skew);
            % indxSKEW=find(skewness(XDupdate')>Th_Skew);
            indxSkewness=makerowvector(indxSKEW);
            % Peaks Ratio of the positive skew signals *********
            % indxPeakRatio=find(ABratio>1); %!!!!!!!!! [ IGNORED ]
            % Make Row Vectors [1xN]:-------------------------
            SNRindx=makerowvector(SNRindx);
            % indxPeakRatio=makerowvector(indxPeakRatio);
            % Get Non-repeated indexes:
            Accepted_index=unique([SNRindx,indxSkewness]);
            Rejected_index=setdiff(1:Ns,Accepted_index);
            if isempty(Rejected_index)
                fprintf('\n\n\n\n > > > > Artifacts Distortion ALERT\n\n\n\n')
            end
            %%% Reject Some False Positives: **************************************
            ActiveNeurons=[];
            if ~isempty(Accepted_index)
                % Get Threshold
                Th_SNR =min(SNRbyWT(Accepted_index));
                % Accepted_index=find(SNRbyWT>=Th_SNR); % UPDATE ACCEPTED
                Th_Skew=get_threshold_pdf(SkewSignal,Accepted_index,Rejected_index);
                % Lil' Fix->
                if Th_Skew<0
                    Th_Skew=0;
                end
                % Apply Threshold
                AcceptedINDX=unique([makerowvector(Accepted_index(SNRbyWT(Accepted_index)>=Th_SNR)),...
                    makerowvector(Accepted_index(SkewSignal(Accepted_index)>Th_Skew))]);
                
                RejectedINDX=setdiff(1:Ns,AcceptedINDX);
                %%% Response Funtion ***********************************************
                [FR(AcceptedINDX,:),~,TAUS(AcceptedINDX,:)]=AR_Estimation(XDupdate(AcceptedINDX,:),p,fs,L,taus_0);
                for k=1:length(AcceptedINDX)
                    if isempty(findpeaks( FR(AcceptedINDX(k),:) ) )
                        FR(AcceptedINDX(k),:)=-FR(AcceptedINDX(k),:);
                        disp('WARNING: Response Function Missestimation')
                    end
                end
                %%% Sparse Deconvolution *******************************************
                [DRIVER(AcceptedINDX,:),LAMBDASS(AcceptedINDX)]=maxlambda_finder(XDupdate(AcceptedINDX,:),FR(AcceptedINDX,:));
                % preDRI=DRIVER;
                % KEEP SNR>0 and High+Skewed Signals
                %% Check Driver - - - - - - - - - - - - - - - - - - - - 
                [Dfix(AcceptedINDX,:),XDfix,Xestfix,LambdasFix,IndexesFix,Features]=...
                    analyze_driver_signal(DRIVER(AcceptedINDX,:),...
                    FR(AcceptedINDX,:),XDupdate(AcceptedINDX,:),Xest(AcceptedINDX,:));
                DRIVER=Dfix;
                if ~isempty(IndexesFix)
                    FixedSignals=AcceptedINDX(IndexesFix);
                    LAMBDASS(FixedSignals,1)=LambdasFix';
                    % XDupdate(FixedSignals,:)=XDfix(FixedSignals,:);
                    XDupdate(FixedSignals,:)=XDfix(IndexesFix,:);
                    Xest(FixedSignals,:)=Xestfix(IndexesFix,:);
                    fprintf('\n\n\n > >   Updated Values   < < \n\n\n')
                    % Process NOT OK#######################################
                    [FRr,~,~]=AR_Estimation(XDupdate(RejectedINDX,:),p,fs,L,taus_0);
                    FR(RejectedINDX,:)=FRr;
                    
                    for k=1:length(RejectedINDX)
                        if isempty(findpeaks( FRr(k,:) ) )
                            FRr(k,:)=-FRr(k,:);
                            disp('WARNING: Response Function Missestimation')
                        end
                    end
                    [DRIVERr,LAMBDASSr]=maxlambda_finder(XDupdate(RejectedINDX,:),FRr,1);
                    LAMBDASS(RejectedINDX)=LAMBDASSr;
                    [DRfix,XDRfix,XestRfix,LambdasRFix,IndexesFixR,FeaturesR]=analyze_driver_signal(DRIVERr,FRr,XDupdate(RejectedINDX,:),Xest(RejectedINDX,:));
                    if ~isempty(IndexesFixR)
                        FixedNOSignals=RejectedINDX(IndexesFixR);
                        LAMBDASS(FixedNOSignals)=LambdasRFix;
                        DRIVER(RejectedINDX,:)=DRfix;
                        XDupdate(RejectedINDX,:)=XDRfix;
                        Xest(RejectedINDX,:)=XestRfix;
                    end
                end
                % Reanalize Sparse Signal to get Active Signals
                disp('Detected Neurons: ')
                for k=1:Ns 
                    xsk=sparse_convolution(DRIVER(k,:),FR(k,:));
                    if ~isempty(xsk(xsk>=std(xsk-XDupdate(k,:)')))
                        ActiveNeurons=[ActiveNeurons;k];
                        
                        fprintf('%i,',ActiveNeurons(end));
                    end
                end
                % Get Active Neurons
                % ActiveNeurons=find( sum(DRIVER,2)~=0 );
            else
                AcceptedINDX=[];
                RejectedINDX=setdiff(1:Ns,AcceptedINDX);
                % DRIVER=[]; % Already Initialized in Zeros
                % FR=[];
                % LAMBDASS=[];
                disp('             *********************' )
                disp('             *********************' )
                disp('             ******PURE NOISE ****' )
                disp('             *********************' )
                disp('             *********************' )
            end
            InactiveNeurons=setdiff(1:Ns,ActiveNeurons);
        end % END WHILE of IndexesFix
        %% Checkin Lambdas
        % if numel(ActiveNeurons)==Ns
        %    smallacceptedlamb=LAMBDASS(ActiveNeurons)<1;
        %    histogram(LAMBDASS(ActiveNeurons(smallacceptedlamb)),numel(find(smallacceptedlamb)))
        %    smallrejectedlamb=LAMBDASS(ActiveNeurons)<1;
        % end
        %% GET RASTER *****************************************************
        % TotalCells=length(XY);
        switch RasterAlgorithm
            case 'Driver'
                Raster=get_raster(1,DRIVER,ActiveNeurons); % DRIVER
            case 'Derivative'
                Raster=get_raster(3,DRIVER,ActiveNeurons,FR); % Derivative
            case 'OOPSI'
                Raster=get_raster(2,DRIVER,ActiveNeurons,FR,fs,Xest); % OOPSI
        end
        % Examples---------------------------------------------------------
        % DRIVER:
        % R1=get_raster(1,Dfix,ActiveNeurons); 
        % OOPSI:
        % R2=get_raster(2,XDupdate(ActiveNeurons,:),ActiveNeurons,TAUS,fs,Xest(ActiveNeurons,:));
        % DERIVATIVE (cleaned up signal)
        % R3=get_raster(3,Dfix,ActiveNeurons,FR);
        % _________________________________________________________________
        % Results Monitor Figure ******************************************
        % figureMonitor=gcf;
        % figureMonitor.Name=[Experiment(2:end),' ',Names_Conditions{i},' vid:',num2str(j)];
        % drawnow;
        
        % Cells to save PREPROCESSING ####################################
        DETSIGNALS{j,i}=XDupdate;       % Detrended Signals         *
        ESTSIGNALS{j,i}=Xest;           % Wavelet Denoised          *
        SNRwavelet{j,i}=SNRbyWT;        % Empirical SNR             
        Responses{j,i}=FR;              % Fluorophore Responses     
        TAUSall{j,i}=TAUS;              % [taus {rise, fall},gain]  
        preDRIVE{j,i}=DRIVER;           % Drivers
        preLAMBDAS{j,i}=LAMBDASS;       % lambda Parameter          
        isSIGNALS{j,i}=ActiveNeurons;   % DETECTTED Signals
        notSIGNALS{j,i}=InactiveNeurons;% UNDETECTED Signals
        RASTER{j,i}=Raster;             % Preliminar Raster         
        
        % Table Data For Processing Log Details ##########################
        TimeProcessing=toc;             % Processing Latency [s]
        T=table( dyename,{num2str(fs)},{num2str(length(ActiveNeurons))},...
            {num2str(Frames)},{num2str(Th_SNR)},{num2str(Th_Skew)},...
            {num2str(TimeProcessing,2)} );

        T.Properties.VariableNames=ColumnNames;
        % Save Table in Resume Tables of the Algorithm Latency*********
        if isdir([FileDirSave,'\Resume Tables'])
            writetable(T,[FileDirSave,'\Resume Tables',[Experiment,'-',Names_Conditions{i}],'.csv'],...
                'Delimiter',',','QuoteStrings',true);
            disp(['Saved Table Resume: ',Experiment,'-',Names_Conditions{i}])
        else % Create Directory
            disp('Directory >Resume Tables< created')
            mkdir([FileDirSave,'\Resume Tables']);
            writetable(T,[FileDirSave,'\Resume Tables',[Experiment,'-',Names_Conditions{i}],'.csv'],...
                'Delimiter',',','QuoteStrings',true);
            disp('Resume Tables Direcotry Created');
            disp(['Saved Table Resume: ',Experiment,'-',Names_Conditions{i}])
        end
        % ARcoeffcients{j,i}=ARc;                 % Autoregressive Coefficients
        % SIGNALSclean{j,i}=X_SPARSE;             % Cleansignals
        % LAMBDASpro{j,i}=LAMBDASproc;            % Sparse Parameter              * To Datasheet
        % SNRs{j,i}=SNRbySD;                      % Signal to Noise Ratio [dB]    * To Datasheet lambdas<1
        % DRIVERs{j,i}=DRIVERSpro;                % Driver Signals
        % QoVid{j,i}=QoV;                         % QUality of Videos  
    end
    disp(' |0|||||||||||||||||||||0||||||||||||||||||||||||0||||||||||||||| ')
    disp(' |||||0||||||||||||||||0||||DATA PROCESSED|||||||||||||||||0||||| ')
    disp(' |||0||||||||||||0||||||||||||||||0|||||||||||||||||||||||||||||| ')
end
%% SAVING(2) Processed Data & Feature Extraction |  Resume Table 
% Save Auto-Processed DATA * * * * * * * * * * * * * * * * * * * * * * * * 
save([FileDirSave,'\Processed Data',Experiment,'.mat'],'DETSIGNALS','ESTSIGNALS',...
    'SNRwavelet','SIGNALSclean','SNRlambda','RasterAlgorithm',...
    'preDRIVE','preLAMBDAS','TAUSall','RASTER','isSIGNALS','notSIGNALS',...
    'Responses','dyename','-append');
disp('Updated: {Feature-Extraction} DATA')
%% Sort & Clean Rasters ***************************************************
% make it nested function--->
% Sort by Activation in each Condition:
[New_Index,Raster_Condition,RASTER_WHOLE]=SortNeuronsCondition(RASTER);
% Plot_Raster_V(RASTER_WHOLE(New_Index,:),fs);
RASTER_WHOLE_Clean=RASTER_WHOLE(New_Index,:);
XY_clean=XY(New_Index,:);
% Clean Raster and Coordinates
ActiveNeurons=find(sum(RASTER_WHOLE_Clean,2)>0);                % INDEX of Active NEURONS only
RASTER_WHOLE_Clean=RASTER_WHOLE_Clean(ActiveNeurons,:);
XY_clean=XY_clean(ActiveNeurons,:);                             % Clean Coordinates
%% PLOT RESULTS
Plot_Raster_V(RASTER_WHOLE_Clean,fs);                           % Clean Whole Raster
set(gcf,'Name',['ID: ',Experiment(2:end),' pre-processing'],'NumberTitle','off')

Label_Condition_Raster(Names_Conditions,Raster_Condition,fs);   % Labels
%% SAVE ReSULTS
save([FileDirSave,'\Processed Data',Experiment,'.mat'],'New_Index','Raster_Condition',...
    'RASTER_WHOLE_Clean','XY_clean','-append');
disp('Saved Sorted Raster Intel')

%% Visual Inpspection & Manual Processing ********************************* GREAT!
% Visual/Manual Proesssing Controler: (-+) & (--)
mf=msgbox({'For Visual Inspection of Detected Ca++ Transients';'Type: ';...
    '>>Detected_Visual_Inspection';' ';'And for Undetected Ca++ Transients';...
    '>>Undetected_Visual_Inspection'});
VisualInspector=[false,false];
waitfor(mf); delete(mf);
%% END OF THE WORLD**************************************************   