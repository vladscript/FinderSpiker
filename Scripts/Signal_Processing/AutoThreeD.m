% Automatic Denoise Detrend & Deconvolution 
% Accepted Signals:
% SNR>0 and Highly +Skewed Signals
function [ESTSIGNALS,SNRwavelet,SNRindx,TAUSall,isSIGNALS,...
    notSIGNALS]=AutoThreeD(NumberofVideos,L,p,taus_0,dyename,FolderNameRT)
%% Call for global variables
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
% Init OutputVariables:
ESTSIGNALS=cell(size(SIGNALS));
SNRwavelet=ESTSIGNALS; TAUSall=ESTSIGNALS;
isSIGNALS=ESTSIGNALS; notSIGNALS=ESTSIGNALS;
%%
NC=size(SIGNALS,2);
% Processing Resume Table Columns:
ColumnNames={'Fluo_Dye','f_s','DetectedCells','Frames',...
        'minSNR','minSkewness','TimeProcessing'};
disp('**************** +++[Processing]+++ *************************')
for i=1:NC
    for j=1:str2double(NumberofVideos{i})
        % Initialize and Read Data
        X=SIGNALS{j,i};
        X_SPARSE=zeros(size(X));
        [Ns,Frames] = size(X);
        FR=zeros(Ns,round(L*fs));
        TAUS=zeros(Ns,3);
        DRIVER=zeros(Ns,Frames);
        LAMBDASS=zeros(Ns,1);
        Dfix=zeros(Ns,Frames);
        % Display Info **********************
        disp(['               [>> Condition: ',num2str(i),'/',num2str(NC),']']);
        disp(['               [>> Video    : ',num2str(j),'/',NumberofVideos{i},']']);
        tic
        % Detrending ******************************************************
        % XD=only_detrending(X);
        % Denoising & Feature Extration ***********************************
        [Xest,SNRbyWT,SkewSignal,~,SkewNoise,XDupdate]=denoise_wavelet(X);
        %% DRIVER ANALYSIS
        IndexesFix=1;       % TO enter to the WhileLoop
        FixedSignals=[];    % 
        FixedNOSignals=[];  % 
        aux=0;              % Loop COunter 
        while and(~isempty(IndexesFix),aux<2)
            aux=aux+1;
            if ~isempty(FixedSignals)
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
            SNRindx=makerowvector(SNRindx);
            % Skew PDFs ****************************************
            % [skew'noise PDF is VERY like RANDN(Ns;Frames)]
            Th_Skew=max(SkewNoise);
            indxSKEW=find(SkewSignal>Th_Skew);
            % indxSKEW=find(skewness(XDupdate')>Th_Skew);
            indxSkewness=makerowvector(indxSKEW);
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
                % KEEP SNR>0 and Highly +Skewed Signals
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
                % Check Sparse Signal to get Active Signals
                disp('Detected Neurons: ')
                for k=1:Ns 
                    xsk=sparse_convolution(DRIVER(k,:),FR(k,:));
                    X_SPARSE(k,:)=xsk';
                    % NOISEX=Xest-XDupdate;
                    if ~isempty(xsk(xsk>=std(xsk-XDupdate(k,:)')))
                        ActiveNeurons=[ActiveNeurons;k];
                        fprintf('%i,',ActiveNeurons(end));
                    end
                end
                % Reprocess to better fit for the Active Ones Only
                [LAMBDASS,X_SPARSE,DRIVER]=fit_sparse(X_SPARSE,Xest,XDupdate,LAMBDASS,FR,DRIVER,ActiveNeurons);
                
            else
                disp('             *********************' )
                disp('             *********************' )
                disp('             ******PURE NOISE ****' )
                disp('             *********************' )
                disp('             *********************' )
            end
            InactiveNeurons=setdiff(1:Ns,ActiveNeurons);
        end % END WHILE of IndexesFix
        %% Checkin Lambdas to FIT BEST
        % % Reprocess to better fit for the Active Ones Only
        [LAMBDASS,X_SPARSE,DRIVER]=fit_sparse(X_SPARSE,Xest,XDupdate,LAMBDASS,FR,DRIVER,ActiveNeurons);
        
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
        
        % Cells to save PREPROCESSING ####################################
        DETSIGNALS{j,i}=XDupdate;       % Detrended Signals
        ESTSIGNALS{j,i}=Xest;           % Wavelet Denoised
        SNRwavelet{j,i}=SNRbyWT;        % Empirical SNR
        Responses{j,i}=FR;              % Fluorophore Responses
        TAUSall{j,i}=TAUS;              % [taus {rise, fall},gain]
        preDRIVE{j,i}=DRIVER;           % Drivers
        preLAMBDAS{j,i}=LAMBDASS;       % lambda Parameter
        SIGNALSclean{j,i}=X_SPARSE;     % Cleansignals
        isSIGNALS{j,i}=ActiveNeurons;   % DETECTTED Signals
        notSIGNALS{j,i}=InactiveNeurons;% UNDETECTED Signals
        RASTER{j,i}=Raster;             % Preliminar Raster
                                        % Signal Noise Ration:
        SNRlambda{j,i}=10*log(var(XDupdate')./var(X_SPARSE')); 
        
        % Table Data For Processing Log Details ##########################
        if ~exist('Th_SNR','var')
            Th_SNR=max(SNRbyWT);
        end
        TimeProcessing=toc;             % Processing Latency [s]
        % Processign Summaries:
        T=table( dyename,{num2str(fs)},{num2str(length(ActiveNeurons))},...
            {num2str(Frames)},{num2str(Th_SNR)},{num2str(Th_Skew)},...
            {num2str(TimeProcessing,2)} );

        T.Properties.VariableNames=ColumnNames;
        % Save Table in Resume Tables of the Algorithm Latency*********
        FileDirSave=getDir2Save();
        if ~isdir([FileDirSave,FolderNameRT])
            mkdir([FileDirSave,FolderNameRT]);
            fprintf('Directory [%s] created',FolderNameRT)
        end
        writetable(T,[FileDirSave,FolderNameRT,[Experiment,'-',Names_Conditions{i}],'.csv'],...
                        'Delimiter',',','QuoteStrings',true);
        disp(['Saved Table Resume: ',Experiment,'-',Names_Conditions{i}])
    end
    disp(' |0|||||||||||||||||||||0||||||||||||||||||||||||0||||||||||||||| ')
    disp(' |||||0||||||||||||||||0||||DATA PROCESSED|||||||||||||||||0||||| ')
    disp(' |||0||||||||||||0||||||||||||||||0|||||||||||||||||||||||||||||| ')
end