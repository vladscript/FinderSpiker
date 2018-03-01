% Function to Manually Detect Calcium Transients
% INPUT
%   isSIGNALS:          Indexes of (detected) Active Neurons
%   SIGNALS             Raw Fluuorescence
%   DETSIGNALS          Detrended FLuorescence
%   preDRIVE            Preliminar Drive
%   preLAMBDAS          Preliminar Lambdas
%   RASTER              Detected Raster ACtivity
%   Responses           Estimated Dye Response
%   Names_Conditions    Names of Conditions
%   fs                  Sampling Frequency
% OUTPUT
%   isSIGNALSOK:        Cell of indexes of detected SIGNALS
%   SIGNALSclean:       Cell of Estimated Signals
%   DRIVEROK:           Cell of driver functions
%   RASTEROK:           Cell of raster neural activity
%   LAMBDASSpro         Cell of Fine Lambdas
%   SNRlambda           Call of SNR from Sparse Analysis
%   OddsMatrix:       Quantification of +P & -N
%           'Frame of Artifacts' maybe
% OPERATION ---------------------------------------------------
% Display Figure and use arrows to :
% Move          Forward  Backward: <-left arrow| right arrow ->
% Deconvolve    More |Less Sparse: up arrow|down arrow
% Click On Raster or Driver Plot to Remove Artifacts
function [isSIGNALSOK,SIGNALSclean,DRIVEROK,RASTEROK,...
        LAMBDASSpro,SNRlambda,OddsMatrix]=...
        Manual_Driver_Raster_Magic(isSIGNALS,SIGNALS,DETSIGNALS,...
        preDRIVE,preLAMBDAS,RASTER,Responses,Names_Conditions,fs)
%% Setup: Initialize Outputs for all Conditions & Videos
isSIGNALSOK=cell(size(isSIGNALS));   % Detected INDEXES
SIGNALSclean=cell(size(DETSIGNALS)); % Detected clean SIGNALS
DRIVEROK=cell(size(preDRIVE));       % clean Drivers
RASTEROK=cell(size(RASTER));         % clean Raster
LAMBDASSpro=cell(size(preLAMBDAS));  % Modified lambdas
SNRlambda=cell(size(preLAMBDAS));    % Sparse Empirical SNRs
delta_lambda=0.25;                   % \lambda Step lambda_next=(1+_deltalambda)*lambda_present
% addpath SpaRSA;                      % Directory of Deconvolution Software
[NVmax,NC]=size(DETSIGNALS);         % NC: N Conditions
OddsMatrix=cell(NVmax,NC);           % Odds Matrix: {++,-+,+-,--}
%% Review : Main Figure: Initialize ++++++++++++++++++++++++++++++++++++++
checksignals=figure('menubar','none','numbertitle','off',...
            'position',[46 42 1220 650],...
            'keypressfcn',@manual_processing_ctrl);

% Define colormap
CM=hot(20);         % Choose Colors
CM=CM(end:-1:1,:);  % Turn up-down
colormap(CM);       % Make Color map
% Initialize AXIS and PLOTS (with randomw stuff) ##########################
ax1=subplot(3,2,1); % Raw Data - - - - - - - - - - - - - -      
plotsignal=plot(randn(1,10),'Parent',ax1);
plotsignal.Parent.YLabel.String='Raw \Delta F/F_{0}';
plotsignal.Parent.XTick=[];
ax2=subplot(3,2,3); % Detrended Data - - - - - - - - - - -
% plotdetrended=plot(randn(1,10),'Parent',ax2);
set(ax2,'NextPlot','replacechildren');
% plotdetrended.Parent.YLabel.String='Detrended';
% plotdetrended.Parent.XTick=[];
ax3=subplot(3,2,5); %  Driver - - - - - - - - - - - - - - - 
plotdriver=bar(randn(1,10),'Parent',ax3);
plotdriver.Parent.YLabel.String='Activity Driver';
set(plotdriver,'ButtonDownFcn',@remove_driver)
ax4=subplot(3,2,[2,4]); % RASTER  - - - - - - - - - - - - - - - 
plotraster=imagesc(randn(10,10),'Parent',ax4);
plotraster.Parent.Title.String='RASTER';
plotraster.Parent.XTick=[];
plotraster.Parent.Box='off';
set(plotraster,'ButtonDownFcn',@remove_columns)
% plotraster.Parent.Box='off';
ax5=subplot(3,2,6); % CoAc Signal  - - - - - - - - - - - - - - -
plotcoac=plot(randn(1,10),'Parent',ax5);
plotcoac.Parent.YLabel.String='CoAc';
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Setup for Processing Noisy Records
L=30;              % seconds  of the Fluorophore Response
p=3;               % AR(p) initial AR order        
taus_0= [.75,2,1]; % starting values for taus
addpath SpaRSA;         % Directory of Deconvolution Software
%% CHECK FALSE POSITIVE ##################################################
for i=1:NC                      % Condition Loop
    for j=1:NVmax               % Video Loop
        checksignals.Name=['False (+) Detection | Condition: ',Names_Conditions{i},' Vid: ',num2str(j)];
        if ~isempty(SIGNALS{j,i})    % Check If There's Video 
            % Read Indexes and Cells:
            [Cells,Frames]=size(SIGNALS{j,i});      % ALL NCells & Frames
            isSIGNAL=isSIGNALS{j,i};                % Detected Signals
            % Check if there're Detected Signals (!)
            if ~isempty(isSIGNAL)    % Check If There's Signals Detected
                preneuron=0;                            % Auxiliar
                % Initialize Outputs (1/2)              Zero Initialization
                SIGNALSclean{j,i}=zeros(Cells,Frames);
                DRIVEROK{j,i}=zeros(Cells,Frames);
                RASTEROK{j,i}=zeros(Cells,Frames);
                LAMBDASSpro{j,i}=zeros(Cells,1);
                SNRlambda{j,i}=zeros(Cells,1);
                % Read Data ++++++++++++++++++++++++++++++++++++++++++
                X=SIGNALS{j,i}(isSIGNAL,:);     % Raw Fluorescence {only-detected}
                FR = Responses{j,i};            % Response Funtions {only-detected}
                XD=DETSIGNALS{j,i}(isSIGNAL,:); % Detrended Fluorescence {only-detected}
                D=preDRIVE{j,i};                % Driver Signals  {only-detected}
                lambdass=preLAMBDAS{j,i};       % List of lambdas {only-detected}
                R=RASTER{j,i};                  % Raster        {ALL CELLS}
                X_SPARSE=zeros(size(X));        % Initialize Sparse Clean Signals
                snrS=zeros(size(lambdass));     % Initialize Output SNRs
                for k=1:length(isSIGNAL)        % {only-detected}
                    xd=XD(k,:);
                    d=D(k,:);
                    r=FR(k,:);
                    x_sparse=sparse_convolution(d,r);
                    X_SPARSE(k,:)=x_sparse;     % Sparse Signals
                    snrc=10*log(var(x_sparse)/var(xd'-x_sparse));
                    snrS(k)=snrc;               % SNR Sparse
                end
                axisdetren=[0,Frames,min(min(XD)),max(max(XD))]; % Axis for detrended Signal
                % Initialize Outputs (2/2)              Automatic Initialization
                SIGNALSclean{j,i}(isSIGNAL,:)=X_SPARSE;
                DRIVEROK{j,i}(isSIGNAL,:)=D;
                RASTEROK{j,i}=R;
                LAMBDASSpro{j,i}(isSIGNAL)=lambdass;
                SNRlambda{j,i}(isSIGNAL)=snrS;
                % Get DATA Initialization +++++++++++++++++++++++++++++++++++
                indx_neuron=1;      % Intialize First Cell
                neuron=isSIGNAL(indx_neuron);                   % Neuron
                x=X(indx_neuron,:);                             % Raw Signal
                r=FR(indx_neuron,:);                            % Response
                xd=XD(indx_neuron,:);                           % Detrended
                d=D(indx_neuron,:);                             % Driver
                daux=d; daux(d<=0)=0;                           % Positve Driver
                x_sparse=X_SPARSE(indx_neuron,:);               % Clean Signal
                snrc=snrS(indx_neuron);                         % SNR Sparse
                lambda=lambdass(indx_neuron);                   % Sparse Paramter
                d=daux;                                         % Driver++
                % Plot Data +++++++++++++++++++++++++++++++++++++++++++++++++
                Plot_Data_Now();        % Where Magic Happens
                % Pause **************************************************
                if indx_neuron<=length(isSIGNAL)
                    % Maybe OK button
                    pause;
                    disp('Another Video in the Wall')
                end
                % Ouput # # # # # # # # # # # # # # # # # # # # # # # # # # #
                % ReviewedSignals
                activeINDX=find(or(sum(R(isSIGNAL,:),2)~=0,snrS>-Inf)); % Indexes of Real Active
                yesisSIGNAL=isSIGNAL( activeINDX );                               % Active Neurons
                isSIGNALSOK{j,i}=yesisSIGNAL;                                     % Save AN
                SIGNALSclean{j,i}(yesisSIGNAL,:)=X_SPARSE(activeINDX,:);          % Update Clean
                DRIVEROK{j,i}(yesisSIGNAL,:)=D(activeINDX,:);                     % Driver Clean
                RASTEROK{j,i}=R(activeINDX,:);                                    % Raster
                LAMBDASSpro{j,i}(yesisSIGNAL)=lambdass(activeINDX);               % Lambdas Fit
                SNRlambda{j,i}(yesisSIGNAL)=snrS(activeINDX);                                  % Sparse SNR
                % Detect False (+) ------------------------------------------
                TruePositive=isSIGNAL(activeINDX);                                % ++
                FalsePositive=isSIGNAL(setdiff(1:length(isSIGNAL),activeINDX));   % -+
                OddsMatrix{j,i}.TruePositive=TruePositive;          % ++
                OddsMatrix{j,i}.FalsePositive=FalsePositive;        % -+
                OddsMatrix{j,i}.FalseNegative=[];                   % -- initialize emtpty
                OddsMatrix{j,i}.TrueNegative=[];                    % +- initialize emtpty
            end
        end
    end
end
close(checksignals);
%% Review : Main Figure: RE-Initialize ++++++++++++++++++++++++++++++++++++++
checksignals=figure('menubar','none','numbertitle','off',...
            'position',[86 42 1220 650],...
            'keypressfcn',@manual_processing_ctrl);
colormap(CM);      % Make Color map
% Initialize AXIS and PLOTS (with randomw stuff) ##########################
ax1=subplot(3,2,1); % Raw Data - - - - - - - - - - - - - -      
plotsignal=plot(randn(1,10),'Parent',ax1);
plotsignal.Parent.YLabel.String='Raw \Delta F/F_{0}';
plotsignal.Parent.XTick=[];
ax2=subplot(3,2,3); % Detrended Data - - - - - - - - - - -
% plotdetrended=plot(randn(1,10),'Parent',ax2);
set(ax2,'NextPlot','replacechildren');
% plotdetrended.Parent.YLabel.String='Detrended';
% plotdetrended.Parent.XTick=[];
ax3=subplot(3,2,5); %  Driver - - - - - - - - - - - - - - - 
plotdriver=bar(randn(1,10),'Parent',ax3);
plotdriver.Parent.YLabel.String='Activity Driver';
set(plotdriver,'ButtonDownFcn',@remove_driver)
ax4=subplot(3,2,[2,4]); % RASTER  - - - - - - - - - - - - - - - 
plotraster=imagesc(randn(10,10),'Parent',ax4);
plotraster.Parent.Title.String='RASTER';
plotraster.Parent.XTick=[];
plotraster.Parent.Box='off';
set(plotraster,'ButtonDownFcn',@remove_columns)
% plotraster.Parent.Box='off';
ax5=subplot(3,2,6); % CoAc Signal  - - - - - - - - - - - - - - -
plotcoac=plot(randn(1,10),'Parent',ax5);
plotcoac.Parent.YLabel.String='CoAc';
%% CHECK FALSE NEGATIVE ###################################################
for i=1:NC                      % Condition Loop
    for j=1:NVmax               % Video Loop
        if ~isempty(SIGNALS{j,i})    % Check If There's Video 
            checksignals.Name=['False (-) Detection | Condition: ',Names_Conditions{i},' Vid: ',num2str(j)];
            % Read Indexes and Cells +++++++++++++++++++++++++++++++++++
            [Cells,Frames]=size(SIGNALS{j,i});      % NCells & Frames
            isSIGNAL=isSIGNALS{j,i};                % Detected Signals
            notSIGNAL=setdiff(1:Cells,isSIGNAL);    % Detected Noisy Records
            preneuron=0;                            % Auxiliar 
            disp('              <[ Process Noisy Detected Records ]>');
            % Read DATA +++++++++++++++++++++++++++++++++++++++++++++++++++
            isSIGNAL=notSIGNAL;              % Posible False (-)
            X=SIGNALS{j,i}(notSIGNAL,:);     % Raw Fluorescence {NOT-detected}
            XD=DETSIGNALS{j,i}(notSIGNAL,:); % Detrended Fluorescence {NOT-detected}
            R=RASTER{j,i};                   % Raster        {ALL CELLS}
            X_SPARSE=zeros(size(X));         % Initialize Sparse Clean Signals
            axisdetren=[0,Frames,min(min(XD)),max(max(XD))]; % Axis for Detrended Signals
            % Process Noise Record +++++++++++++++++++++++++++++++++++++++
            [FR,~,~]=AR_Estimation(XD,p,fs,L,taus_0);
            for k=1:length(notSIGNAL)
                if isempty(findpeaks( FR(k,:) ) )
                    FR(k,:)=-FR(k,:);
                    disp('WARNING: Response Function Missestimation')
                end
            end
            [D,lambdass]=maxlambda_finder(XD,FR);
            snrS=zeros(size(lambdass));      % Initialize Output SNRs
            % Clean Drivers (only +++Drivers) ++++++++++++++++++++++++++++
            ThDriver=abs(min(D'));
            [Nd,~]=size(D);
            for n=1:Nd
                D(n,D(n,:)<ThDriver(n))=0;
            end
            % Get Sparse Signals and SNRs
            for k=1:length(isSIGNAL)        % {only-detected}
                xd=XD(k,:);
                d=D(k,:);
                r=FR(k,:);
                x_sparse=sparse_convolution(d,r);
                X_SPARSE(k,:)=x_sparse;     % Sparse Signals
                snrc=10*log(var(x_sparse)/var(xd'-x_sparse));   
                snrS(k)=snrc;               % SNR Sparse
            end
            % Get Raster +++++++++++++++++++++++++++++++++++++++++++++++++
            for n=1:length(notSIGNAL)
                [~,Np]=findpeaks(D(n,:)); % way too clean
                % [~,Np]=find(D(n,:)>0);  % alright
                R(notSIGNAL(n),Np)=1;
            end
            % Initialize Outputs               Automatic Initialization
            SIGNALSclean{j,i}(isSIGNAL,:)=X_SPARSE;
            DRIVEROK{j,i}(isSIGNAL,:)=D;
            RASTEROK{j,i}=R;
            LAMBDASSpro{j,i}(isSIGNAL)=lambdass;
            SNRlambda{j,i}(isSIGNAL)=snrS;
             % Get DATA Initialization +++++++++++++++++++++++++++++++++++
            indx_neuron=1;      % Intialize First Cell
            if ~isempty(isSIGNAL)    % Check If There were Non-Detected Signals 
                neuron=isSIGNAL(indx_neuron);                   % Neuron
                x=X(indx_neuron,:);                             % Raw Signal
                r=FR(indx_neuron,:);                            % Response
                xd=XD(indx_neuron,:);                           % Detrended
                d=D(indx_neuron,:);                             % Driver
                daux=d; daux(d<=0)=0;                           % Positve Driver
                x_sparse=X_SPARSE(indx_neuron,:);               % Clean Signal
                snrc=snrS(indx_neuron);                         % SNR
                lambda=lambdass(indx_neuron);  
                d=daux;
                % Plot Data +++++++++++++++++++++++++++++++++++++++++++++++++
                Plot_Data_Now();        % Where Magic Happens
                % Pause **************************************************
                if indx_neuron<=length(isSIGNAL)
                    % Maybe OK button
                    pause;
                    disp('Another Video in the Wall')
                end
                % Ouput # # # # # # # # # # # # # # # # # # # # # # # # # # #
                % Reviewed Signals
                recoveredINDX=find( or(sum(R(isSIGNAL,:),2)~=0,snrS>-Inf) );
                recoSIGNAL=isSIGNAL( recoveredINDX );
                isSIGNALSOK{j,i}=union(isSIGNALSOK{j,i},recoSIGNAL);
                SIGNALSclean{j,i}(recoSIGNAL,:)=X_SPARSE(recoveredINDX,:);
                DRIVEROK{j,i}(recoSIGNAL,:)=D(recoveredINDX,:);
                RASTEROK{j,i}=R;
                LAMBDASSpro{j,i}(recoSIGNAL)=lambdass(recoveredINDX);
                SNRlambda{j,i}(recoSIGNAL)=snrS(recoveredINDX);
                % Detect False (+) ------------------------------------------
                FalseNegative=isSIGNAL(recoveredINDX);                               % -- 
                TrueNegative=isSIGNAL(setdiff(1:length(isSIGNAL),recoveredINDX));   % +-
                OddsMatrix{j,i}.FalseNegative=FalseNegative;        % --
                OddsMatrix{j,i}.TrueNegative=TrueNegative;          % +-
                if isempty(isSIGNALS{j,i})
                    OddsMatrix{j,i}.TruePositive=[];          % ++
                    OddsMatrix{j,i}.FalsePositive=[];        % -+
                end
            end
        end
    end
end
close(checksignals);
% update main function output


%% Nested Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% function Set_Name_Axis()
%         % Figure's Name
%         if CHKtype==1
%             checksignals.Name=['False (+) Detection | Condition: ',Names_Conditions{i},' Vid: ',num2str(j)];
%         else
%             checksignals.Name=['False (-) Detection | Condition: ',Names_Conditions{i},' Vid: ',num2str(j)];
%         end
%         % Initialize AXIS;
%         ax1=subplot(3,2,1); % Detrended Data - - - - - - - - - - -      
%         ax1.YLabel.String='Raw \Delta F/F_{0}';
%         ax1.XTick=[];
%         ax2=subplot(3,2,3); % Detrended Data - - - - - - - - - - -
%         ax2.YLabel.String='Detrended \Delta F/F_{0}';
%         ax2.XTick=[];
%         ax3=subplot(3,2,5); %  Driver - - - - - - - - - - - - - - - 
%         ax3.YLabel.String='driver';
%         % ax3.YTick=zeros(1,0);
%         ax4=subplot(3,2,[2,4]); % RASTER  - - - - - - - - - - - - - - - 
%         ax4.Title.String='RASTER';
%         ax4.XTick=[];
%         ax4.Box='off';
%         ax5=subplot(3,2,6); % CoAc Signal  - - - - - - - - - - - - - - -
%         ax5.YLabel.String='CoAc';
% end

function Get_Data()
    neuron=isSIGNAL(indx_neuron);                   % Neuron
    x=X(indx_neuron,:);                             % Raw Signal
    r=FR(indx_neuron,:);                            % Response
    xd=XD(indx_neuron,:);                           % Detrended
    d=D(indx_neuron,:);                             % Driver
    daux=d; daux(d<=0)=0;                           % Positve Driver
    x_sparse=X_SPARSE(indx_neuron,:);               % Clean Signal
    % x_sparse=sparse_convolution(d,r);               % Clean Signal
    snrc=10*log(var(x_sparse)/var(xd-x_sparse));   % SNR
    lambda=lambdass(indx_neuron);  
end


function Plot_Data_Now()
% Funtion to plot data for each VIDEO
        % Raw Data - - - - - - - - - - - - - -  
        plotsignal.YData=x;
        plotsignal.Parent.XLim=[1,Frames];
        plotsignal.Parent.YLim=[min(x),max(x)];
        plotsignal.Parent.XGrid='on';
        plotsignal.Parent.YGrid='on';
        title(plotsignal.Parent,['Neuron: ',num2str(neuron),...
            ' Progress: ',num2str(100*indx_neuron/length(isSIGNAL),3),'%',],...
              'FontSize',8,'FontWeight','normal');
        
        % Detrended Data - - - - - - - - - - -
        % plotdetrended.YData=xd;
        cla(ax2);
        hold(ax2,'on')
        plot(ax2,xd,'b','LineWidth',1)
        plot(ax2,x_sparse,'r','LineWidth',2)
        grid(ax2,'on');
        hold(ax2,'off')
        axis(ax2,axisdetren);
        title(ax2,['SNR=',num2str(snrc),'[dB]'],...
              'FontSize',8,'FontWeight','normal');
        
        %  Driver - - - - - - - - - - - - - - - 
        daux=daux/max(daux);    % Normalize
        plotdriver.YData=daux;
%         if ~isempty(daux(daux>0))||~isinf(snrc)
%             
%             bar(plotdriver.Parent,daux,'FaceColor',[0,1,0.5]);
%         else
%             plotdriver.YData=daux;
%             bar(plotdriver.Parent,daux,'FaceColor','r');
%         end
        axis(plotdriver.Parent,'tight');
        grid(plotdriver.Parent,'on');
        title(plotdriver.Parent,['\lambda=',num2str(lambda)],...
             'FontSize',8,'FontWeight','normal');
        
        % RASTER  - - - - - - - - - - - - - - - 
        % subplot(3,2,[2,4]); % RASTER
        Raux=R*2;
        if ~isempty(Raux(neuron,:)>0)
            Raux(neuron,Raux(neuron,:)>0)=1;
        end
        plotraster.CData=Raux;
        plotraster.Parent.YLim=[0,Cells];
        plotraster.Parent.XLim=[1,Frames];
        % imagesc(plotraster.Parent,Raux);
        % axis(plotraster.Parent,[1,Frames,0,Cells])
        plotraster.Parent.XTick=[];
        plotraster.Parent.YTick=[];
        
        % CoAc Signal  - - - - - - - - - - - - - - -
        plotcoac.YData=sum(R);
        % plot(plotcoac.Parent,sum(R),'LineWidth',2);
        axis(plotcoac.Parent,'tight');
        grid(plotcoac.Parent,'on');
        drawnow;
end

function manual_processing_ctrl(checksignals,~,~)
        key=get(checksignals,'CurrentKey');
        if strcmp(key,'rightarrow')
            disp('Next ------->')
            if indx_neuron<length(isSIGNAL)
                preneuron=isSIGNAL(indx_neuron);
                indx_neuron=indx_neuron+1;
            else
                preneuron=isSIGNAL(indx_neuron);
                indx_neuron=length(isSIGNAL);
            end
%             [neuron,x,xd,daux,r,x_sparse,snrc,lambda]=Get_Data(indx_neuron,...
%                 X,FR,XD,D,lambdass,isSIGNAL); 
            Get_Data;  d=daux;
            disp(neuron)
            Plot_Data_Now;
            % Set_Name_Axis;
            update_data;
            % SHOULD SAVE DATA
        elseif strcmp(key,'leftarrow')
            disp('<-- Previous')
            if indx_neuron==1
                preneuron=isSIGNAL(indx_neuron);
                indx_neuron=1;
            else
                preneuron=isSIGNAL(indx_neuron);
                indx_neuron=indx_neuron-1;
            end
%              [neuron,x,xd,daux,r,x_sparse,snrc,lambda]=Get_Data(indx_neuron,...
%                 X,FR,XD,D,lambdass,isSIGNAL); d=daux;
            Get_Data;  d=daux;
            disp(neuron)
            Plot_Data_Now;
            % Set_Name_Axis
            update_data;
        elseif strcmp(key,'uparrow')
            disp('Increase lambda /\')
            disp('............... |')
            Get_Data;  d=daux;
            lambda=(1+delta_lambda)*lambda;
            [d,x_sparse,~]=magic_sparse_deconvolution(xd,r,lambda);
            daux=d; daux(d<=0)=0;
            R(neuron,:)=0;
            R(neuron,daux>0)=1;
            sparsenoise=xd-x_sparse';
            snrc=10*log(var(x_sparse)/var(sparsenoise)); 
            Plot_Data_Now;
            update_data;
        elseif strcmp(key,'downarrow')
            disp('...............  |')
            disp('Decrease lambda \/')
            Get_Data;  d=daux;
            lambda=(1-delta_lambda)*lambda;
            [d,x_sparse,~]=magic_sparse_deconvolution(xd,r,lambda);
            daux=d; daux(d<=0)=0;
            R(neuron,:)=0;
            R(neuron,daux>0)=1;
            sparsenoise=xd-x_sparse';
            snrc=10*log(var(x_sparse)/var(sparsenoise)); 
            Plot_Data_Now;
            update_data;
        end
 end
 
 function update_data()
        XD(indx_neuron,:)=xd;          
        D(indx_neuron,:)=daux;
        snrS(indx_neuron)=snrc;
        X_SPARSE(indx_neuron,:)=x_sparse;
        lambdass(indx_neuron)=lambda;
        % Missing Raster?
 end
 function remove_driver(~,~)
     % rect: [xmin ymin width height]
    rect=getrect; rect=round(rect);
    frames2del=[rect(1),rect(1)+rect(3)];
    % Limit Window to Delete
    if frames2del(1)<1
        frames2del(1)=1;
    end
    if frames2del(2)>Frames
        frames2del(2)=Frames;
    end
    daux(frames2del(1):frames2del(2))=0;
    R(neuron,frames2del(1):frames2del(2))=0;
    Plot_Data_Now;      % update plot
    update_data;        % saves daux
    sprintf('Deleted Driver Frames: %d : %d',frames2del)
 end

 function remove_columns(~,~)
     % rect: [xmin ymin width height]
    rect=getrect; rect=round(rect);
    frames2del=[rect(1),rect(1)+rect(3)];
    % Limit Window to Delete
    if frames2del(1)<1
        frames2del(1)=1;
    end
    if frames2del(2)>Frames
        frames2del(2)=Frames;
    end
    % daux(frames2del(1):frames2del(2))=0;
    R(:,frames2del(1):frames2del(2))=0;
    Plot_Data_Now;      % update plot
    update_data;        % saves daux
    sprintf('Deleted Raster Frames: %d : %d',frames2del)
 end

end % END OF THE WORLD

%% Nested Funnctions
% Manual Controller Processing ****************
% Use keyboard input for processing signals
% Input
% Current Figure
% Up    Arrow
% Down  Arrow
% Left  Arrow
% Right Arrow
% Output
