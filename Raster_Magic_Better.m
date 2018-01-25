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
%   ExperimentID.mat @Results Folder        Useful Workspace
%   Settings_Automatic_Processing_Log.xls   Features from automatic mode
%   Features_Processgin_Table.xls           Features from manual mode

%   UPDATES
%   ManualMode&SavingFeatures       15/11/17
%   Huge One @Algorithmie           02/10/17 (no se olvida)
%   Using GitHub                    19/01/2017
%% Global Setup ***********************************************************
clc
clear;
close all;

%% ADDING ALLSCRIPTS
ActualDir=pwd;
addpath(genpath([ActualDir,'\Scripts']))

%% Set Default Directory of Experiments
DefaultPath='C:\Users\Vladimir\Documents\Doctorado\Experimentos\';
if exist(DefaultPath,'dir')==0
    DefaultPath=pwd;
end
% Read Sampling Frequency
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
NV=max(str2num(cell2mat(NumberofVideos))); % Max N of Videos
[~,NC]=size(FN);                           % N Conditions
%% Initalization Data Output
SIGNALS=cell(NV,NC);
DETSIGNALS=cell(NV,NC);
ESTSIGNALS=cell(NV,NC);
SNRwavelet=cell(NV,NC);
Responses=cell(NV,NC);
preDRIVE=cell(NV,NC);
preLAMBDAS=cell(NV,NC);
TAUSall=cell(NV,NC);
SIGNALSclean=cell(NV,NC);
RASTER=cell(NV,NC);
isSIGNALS=cell(NV,NC);
LAMBDASpro=cell(NV,NC);
SNRs=cell(NV,NC);
DRIVERs=cell(NV,NC);

%% Load Data *********************************************************


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
%% Save(1) RAW Data * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
% Experiment: Name of the folder
FileDirSave=pwd;
slashes=find(PathName=='\');
Experiment=PathName(slashes(end-1):slashes(end)-1);
save([FileDirSave,'\Results',Experiment,'.mat'],'Experiment','SIGNALS',...
    'Names_Conditions','NumberofVideos','XY','fs','r');

%% SETUP PROCESSING PARAMTERS ******************************
% Setup for Auto Regressive Process Estimation
L=30;              % seconds  of the Fluorophore Response
p=3;               % AR(p) initial AR order        
taus_0= [.75,2,1]; % starting values for taus
% addpath SpaRSA;         % Directory of Deconvolution Software

%% DETRENDING AND SPARSE DECONVOLVE DETECTION
[NV,NC]=size(SIGNALS);
disp('**************** +++[Processing]+++ *************************')
for i=1:NC
% for i=2:NC
    auxj=0;
    for j=1:str2double(NumberofVideos{i})
%     for j=2:str2double(NumberofVideos{i})
        X=SIGNALS{j,i};
        % Display Info **********************
        disp(['               [>> Condition: ',num2str(i),'/',num2str(NC),']']);
        disp(['               [>> Video    : ',num2str(j),'/',NumberofVideos{i},']']);
        tic
        % Detrending ******************************************************
        [Ns,Frames] = size(X);
        XD=only_detrending(X);
        %%% Denoising ******************************************************
        % Main Features to Detect Signal or Noise:
        % [Xest,SNRbyWT,SkewSignal,ABratio,SkewNoise]=denoise_wavelet(XD,X);
        [Xest,SNRbyWT,SkewSignal,ABratio,SkewNoise,XDupdate]=denoise_wavelet(XD);
        %%% Decision Features
        % SNR >0 *******************************************
        SNRindx=find(SNRbyWT>0);
        % Skew PDFs ****************************************
        % [skew'noise PDF is VERY like RANDN(Ns;Frames)]
        Th_Skew=max(SkewNoise);
        indxSKEW=find(SkewSignal>Th_Skew);
        indxSkewness=makerowvector(indxSKEW);
        % Peaks Ratio of the positive skew signals *********
        indxPeakRatio=find(ABratio(indxSKEW)>1); %!!!!!!!!! REMOVED
        % Make Row Vectors [1xN]:-------------------------
        SNRindx=makerowvector(SNRindx);
        indxPeakRatio=makerowvector(indxPeakRatio);
        % Get Non-repeated indexes:
        Accepted_index=unique([SNRindx,indxSkewness]);
        Rejected_index=setdiff(1:Ns,Accepted_index);
        %%% Reject Som False Positives: **************************************
        if ~isempty(Accepted_index)
            % Get Threshold
            Th_SNR =get_threshold_pdf(SNRbyWT,Accepted_index,Rejected_index);
            Accepted_index=find(SNRbyWT>=Th_SNR); % UPDATE ACCEPTED
            Th_Skew=get_threshold_pdf(SkewSignal,Accepted_index,Rejected_index);
            % Apply Threshold
            A=Accepted_index(SNRbyWT(Accepted_index)>Th_SNR);
            B=Accepted_index(SkewSignal(Accepted_index)>Th_Skew);
            A=makerowvector(A);
            B=makerowvector(B);
            AcceptedINDX=unique([A,B]);
            RejectedINDX=setdiff(1:Ns,AcceptedINDX);
            %%% Plot Some PDFs to monitor errors (!)
%             figure; 
%             subplot(211)
%             if and(~isempty(AcceptedINDX),~isempty(RejectedINDX))
%                 [pA,binA]=ksdensity(SNRbyWT(AcceptedINDX),linspace(min(SNRbyWT(AcceptedINDX)),max(SNRbyWT(AcceptedINDX)),100));
%                 [pR,binR]=ksdensity(SNRbyWT(RejectedINDX),linspace(min(SNRbyWT(RejectedINDX)),max(SNRbyWT(RejectedINDX)),100));
%                 plot(binR,pR,'.r','LineWidth',2)
%                 hold on;
%                 plot(binA,pA,'*b','LineWidth',2)
%                 hold off;
%                 legend('Rejected','Accepted','Location','northwest')
%             else
%                 [pO,binO]=ksdensity(SNRbyWT,linspace(min(SNRbyWT),max(SNRbyWT),100));
%                 plot(binO,pO,'.m','LineWidth',1)
%             end
%             axis tight; grid on;
%             title('SNR pdf')
%             subplot(212)
%             if and(~isempty(AcceptedINDX),~isempty(RejectedINDX))
%                 [pA,binA]=ksdensity(SkewSignal(AcceptedINDX),linspace(min(SkewSignal(AcceptedINDX)),max(SkewSignal(AcceptedINDX)),100));
%                 [pR,binR]=ksdensity(SkewSignal(RejectedINDX),linspace(min(SkewSignal(RejectedINDX)),max(SkewSignal(RejectedINDX)),100));
%                 plot(binR,pR,'.r','LineWidth',2)
%                 hold on;
%                 plot(binA,pA,'*b','LineWidth',2)
%                 hold off;
%                 legend('Rejected','Accepted','Location','northwest')
%             else
%                 [pO,binO]=ksdensity(SkewSignal,linspace(min(SkewSignal),max(SkewSignal),100));
%                 plot(binO,pO,'.m','LineWidth',1)
%             end
%             axis tight; grid on;
%             title('Skewness pdf')
            
            %%% Response Funtion ***********************************************
            [FR,~,TAUS]=AR_Estimation(XD(AcceptedINDX,:),p,fs,L,taus_0);
            for k=1:length(AcceptedINDX)
                if isempty(findpeaks( FR(k,:) ) )
                    FR(k,:)=-FR(k,:);
                    disp('WARNING: Response Function Missestimation')
                end
            end
            % Search in MODES of TAUS  -----------------------------------
            % [Rcanon,TAU_PEAK]=canonical_response(TAUS,fs,1);
            % % Make positive Responses {if needed}
            %if isempty(findpeaks(Rcanon))
            %    Rcanon=-Rcanon;
            %    disp('WARNING: Response Function Missestimation')
            %end
            % FR=repmat(Rcanon,numel(AcceptedINDX),1);
            
            %%% Sparse Deconvolution *******************************************
            [DRIVER,LAMBDASS]=maxlambda_finder(XD(AcceptedINDX,:),FR);
            % Ignore Drivers with bigger Negative Drivers than positive ones
            Dindex=find( abs(min(DRIVER'))<abs(max(DRIVER')) );
            % Update Variables:
            D=DRIVER(Dindex,:);
            FR=FR(Dindex,:);
            LAMBDASS=LAMBDASS(Dindex);
            ActiveNeurons=AcceptedINDX(Dindex);
            InactiveNeurons=setdiff(1:Ns,ActiveNeurons);
            % Negative Threshold to clean Drivers ****************************
            ThDriver=abs(min(D'));
            % Clean Drivers (only +++Drivers)
            [Nd,~]=size(D);
            for n=1:Nd
                D(n,D(n,:)<ThDriver(n))=0;
            end
            % Check if zero drivers & Update DATA
            ActiveNeurons=ActiveNeurons( sum(D,2)~=0 );
            InactiveNeurons=setdiff(1:Ns,ActiveNeurons);
            LAMBDASS=LAMBDASS( sum(D,2)~=0 );
            D=D(sum(D,2)~=0, :);
            FR=FR(sum(D,2)~=0,:);
            
        else
            AcceptedINDX=[];
            RejectedINDX=setdiff(1:Ns,AcceptedINDX);
            D=[];
            ActiveNeurons=[];
            FR=[];
            LAMBDASS=[];
            %% plot to monitor results
%             figure; 
%             subplot(211)
%             [pO,binO]=ksdensity(SNRbyWT,linspace(min(SNRbyWT),max(SNRbyWT),100));
%             plot(binO,pO,'.k','LineWidth',1)
%             axis tight; grid on;
%             title('SNR pdf')
%             subplot(212)
%             [pO,binO]=ksdensity(SkewSignal,linspace(min(SkewSignal),max(SkewSignal),100));
%             plot(binO,pO,'.k','LineWidth',1)
%             axis tight; grid on;
%             title('Skewness pdf')
            disp('             *********************' )
            disp('             *********************' )
            disp('             ******PURE NOISE ****' )
            disp('             *********************' )
            disp('             *********************' )
        end
        %%% GET RASTER *****************************************************
        Raster=zeros(size(XD));
        if ~isempty(ActiveNeurons)
            [AN,~]=size(D);
            for n=1:AN
                % [~,Np]=findpeaks(D(n,:)); % way too clean
                [~,Np]=find(D(n,:)>0);      % alright
                Raster(ActiveNeurons(n),Np)=1;
            end
        end
        % Monitor Figure **************************************************
%         figureMonitor=gcf;
%         figureMonitor.Name=[Experiment(2:end),' ',Names_Conditions{i},' vid:',num2str(j)];
%         drawnow;

        % Cells to save PREPROCESSING ####################################
        
        DETSIGNALS{j,i}=XD;             % Detrended Signals
        ESTSIGNALS{j,i}=Xest;           % Wavelet Denoised
        SNRwavelet{j,i}=SNRbyWT;        % Empirical SNR             * To Datasheet
        Responses{j,i}=FR;              % Fluorophore Responses
        TAUSall{j,i}=TAUS;            % [taus {rise, fall},gain]    * To Datasheet
        preDRIVE{j,i}=D;                % Preliminar Drives (+ & -)
        preLAMBDAS{j,i}=LAMBDASS;       % lambda Parameter          * To Datasheet
        isSIGNALS{j,i}=ActiveNeurons;   % Signal Indicator          * To Datasheet       
        RASTER{j,i}=Raster;             % Preliminar Raster
        % Table Data For Processing Log Details ##########################
        TimeProcessing=toc;             % Processing Latency [s]
        T=table( dyename,{Experiment},{num2str(fs)},{Names_Conditions{i}},...
            {num2str(length(ActiveNeurons))},{num2str(Frames)},...
            {num2str(Th_SNR)},{num2str(Th_Skew)}, {num2str(TimeProcessing,2)} );
%         SavetoExcel(T);
        % ARcoeffcients{j,i}=ARc;       % Autoregressive Coefficients
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
 %% SAVING(2)
% Save Auto-Processed DATA * * * * * * * * * * * * * * * * * * * * * * * * 
save([FileDirSave,'\Results',Experiment,'.mat'],'DETSIGNALS','ESTSIGNALS','SNRwavelet',...
    'preDRIVE','preLAMBDAS','TAUSall','RASTER','isSIGNALS','Responses','dyename','-append');

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
set(gcf,'Name',['ID: ',Experiment(2:end)],'NumberTitle','off')

Label_Condition_Raster(Names_Conditions,Raster_Condition,fs);   % Labels
%% SAVE ReSULTS
save([FileDirSave,'\Results',Experiment,'.mat'],'New_Index','Raster_Condition',...
    'RASTER_WHOLE_Clean','XY_clean','-append');
disp('Saved at Results Folder')
% useful:
% addpath SpaRSA; % Directory of Deconvolution Software
% addpath NeuralNetworks; % NeuralNetworks: Ensamble Analysis
%% Visual Inpspection & Manual Processing ********************************* GREAT!
% Ask if so
button = questdlg('Results Inspection?');
if strcmp('Yes',button)
    [NV,NC]=size(isSIGNALS);
    preisSIGNALS=isSIGNALS;
    [isSIGNALSOK,SIGNALSclean,DRIVEROK,RASTEROK,...
    LAMBDASSpro,SNRlambda,OddsMatrix]=Manual_Driver_Raster_Magic(isSIGNALS,SIGNALS,...
    DETSIGNALS,preDRIVE,preLAMBDAS,RASTER,Responses,Names_Conditions,fs);
    
    % Re-Sort WHOLE RASTER 
    [New_Index,Raster_Condition,RASTER_WHOLE]=SortNeuronsCondition(RASTEROK);
    RASTER_WHOLE_Clean=RASTER_WHOLE(New_Index,:);
    XY_clean=XY(New_Index,:);
    % Clean Raster and Coordinates
    TotalActiveNeurons=find(sum(RASTER_WHOLE_Clean,2)>0);                % INDEX of Active NEURONS
    QoE=round(100*length(TotalActiveNeurons)/length(XY),2);
    sprintf('Actual Active Neurons: %d %%',round(100*QoE));
    % Whole Raster
    RASTER_WHOLE_Clean=RASTER_WHOLE_Clean(TotalActiveNeurons,:);
    XY_clean=XY_clean(TotalActiveNeurons,:);                        % Clean Coordinates
    % SEE RESULTS ################################################
    Plot_Raster_V(RASTER_WHOLE_Clean,fs);                           % Clean Whole Raster
    set(gcf,'Name',['ID: ',Experiment(2:end)],'NumberTitle','off')
    Label_Condition_Raster(Names_Conditions,Raster_Condition,fs);   % Labels    
    % Update Results                [ok] ----------------------------------
    % Ask for Directory to save & MAT file to update
    checkname=1;
    while checkname==1
        DefaultPath='C:\Users\Vladimir\Documents\Doctorado\Software\GetTransitum\Calcium Imaging Signal Processing\Results';
            if exist(DefaultPath,'dir')==0
                DefaultPath=pwd; % Current Diretory of MATLAB
            end

        [FileName,PathName] = uigetfile('*.mat',[' Pick the Analysis File ',Experiment],...
            'MultiSelect', 'off',DefaultPath);
        dotindex=find(FileName=='.');
        if strcmp(FileName(1:dotindex-1),Experiment(2:end))
            checkname=0;
            % SAVE DATA
            save([PathName,FileName],'isSIGNALSOK','SIGNALSclean','DRIVEROK','SNRlambda','LAMBDASSpro',...
                'New_Index','Raster_Condition','RASTER_WHOLE_Clean','XY_clean','RASTEROK','-append');
            disp([Experiment,'   -> UPDATED'])

        else
            disp('Not the same Experiment!')
            disp('Try again!')
        end
     end    
    %% Save Processing-Features of RASTER Activity: p(++), p(-+), p(--),p(+-)
    % Get N videos and N conditions
    [Nv,Nc]=size(ESTSIGNALS);
    NNeurons=length(XY);
    % SETUP ACCUMULATIVE ARRAYS
    EXPNAME=[];
    FLUOROPHORE=[];
    CONDNAME={};
    COORD=[];
    RAD=[];
    SNRW=[];
    SNRL=[];
    SKEWDEN=[];
    LAMB=[];
    ODDS=[];
    for c=1:Nc
        for v=1:Nv
            if ~isempty(ESTSIGNALS{v,c})
                EXPNAME=[EXPNAME;repmat(Experiment,NNeurons,1)];
                FLUOROPHORE=[FLUOROPHORE;repmat(dyename{1},NNeurons,1)];
                % CONDNAME=[CONDNAME;repmat(Names_Conditions{c},NNeurons,1)];
                LWord=length(Names_Conditions{c});
                CONDNAME=[CONDNAME;mat2cell(repmat(Names_Conditions{c},NNeurons,1),...
                    ones(NNeurons,1),LWord)];
                COORD=[COORD;XY];
                RAD=[RAD;r];
                SNRW=[SNRW;SNRwavelet{v,c}];
                SNRL=[SNRL;makerowvector(SNRlambda{v,c})'];
                Xden=ESTSIGNALS{v,c}; 
                SKEWDEN=[SKEWDEN;skewness(Xden')']; % Tranpose it to make it row vector
                LAMB=[LAMB;makerowvector(LAMBDASSpro{v,c})']; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ManualMode
                % Odds Matrix Treatment
                Detection=categorical;
                OM=OddsMatrix{v,c};
                Detection(OM.TruePositive)='++';
                Detection(OM.TrueNegative)='+-';
                Detection(OM.FalsePositive)='-+';
                Detection(OM.FalseNegative)='--';
                ODDS=[ODDS;Detection'];     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ManualMode
            end
        end
    end
    Tfeat=table(EXPNAME,FLUOROPHORE,CONDNAME,COORD,RAD,SNRW,SNRL,SKEWDEN,LAMB,ODDS);
    Tfeat.Properties.VariableNames={'Experiment','Dye','Condition','ROIcoordinates',...
        'ROIradius','SNRwavelet','SNRdeconv','SignalSkewness','lambda','Detection'};
    disp('Saving...')
    summary(Tfeat)
    FileProcessingFeatures='Features_Processgin_Table.xls';
    if exist(FileProcessingFeatures,'file')>0 % If the File Already Exists
        disp('Writing on Table ...')
        Tprevious=readtable(FileProcessingFeatures);
        [ColumnsUsed,~]=size(Tprevious);
        corner_indx=['A',num2str(ColumnsUsed+2)];
        writetable(Tfeat,FileProcessingFeatures,'Sheet',1,'Range',corner_indx,...
            'WriteVariableNames',0);
    else                            % New File
        disp('Creating Table Processing Settings...')
        writetable(Tfeat,FileProcessingFeatures,'Sheet',1);
    end
    disp('saved.')
end

%% END OF THE WORLD**************************************************