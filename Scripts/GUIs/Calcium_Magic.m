%% Function to Clean Manuallly Calcium Transients
% Input
% indxSIGNAL: Cell of Indexed Signal as Detected Ca^2+ Transient
% Global Inputs:
%   SIGNALS             Raw Fluuorescence
%   DETSIGNALS          Detrended Fluorescence
%   preDRIVE            Preliminar Drive
%   preLAMBDAS          Preliminar Lambdas
%   RASTER              Detected Raster ACtivity
%   Responses           Estimated Dye Response
%   Names_Conditions    Names of Conditions
% Outputs
%
%   indxSIGNALSOK:      Acutal Detected Signals
% Global OutputModifiied
%   SIGNALSclean:       Sparse Cleaned up Signals
%   SNRs                Signal Noise Ratio
%   preDRIVE            Preliminar Drive
%   preLAMBDAS          Preliminar Lambdas
%   RASTER              Detected Raster ACtivity
%   Responses           Estimated Dye Response

function indxSIGNALSOK = Calcium_Magic(indxSIGNALS)
%% Setup
disp('Loading...')
% Call Global Variables
% global isSIGNALS;
% global notSIGNALS;
global SIGNALS;             % Raw Fluorescence Signals
global DETSIGNALS;          % Detrended Fluorescence Signals
global preDRIVE;            % Automatic Founded Driver Function
global preLAMBDAS;          % Automatic Lambdas Founded
global RASTER;              % Raster Activity
global Responses;           % Fluorescence Response Signals
global RasterAlgorithm;     % Algorithm to Get Raster: Driver,OOPSI, Derivative
global Names_Conditions;    % Names of The Conditions
global SIGNALSclean;        % Sparse Clean Signals
global SNRlambda;           % Signal Noise Ratio by Sparse Deconvolution
global Experiment;          % Experiment ID
Experiment=Experiment(Experiment~='\');
global fs;
global indxSIGNALSOK;       % Indexes of the Detected Signals
% indxSIGNALSOK=cell(size(indxSIGNALS));
% Setup Variables
% Load_Default_Values_SP;
delta_lambda=0.25;          % \lambda Step lambda_next=(1+_deltalambda)*lambda_present
[NVmax,NC]=size(SIGNALS);   % NC: N Conditions
% Get Useful Pair of Indexes i,j
IndxPair=[];
for ihat=1:NC
    for jhat=1:NVmax
       if and(~isempty(SIGNALS{jhat,ihat}),~isempty(indxSIGNALS{jhat,ihat}))
       % if ~isempty(indxSIGNALS{jhat,ihat})
           IndxPair=[IndxPair;jhat,ihat];
       end
    end
end
Npairs=size(IndxPair,1);
% Load_Default_Directories;
%% Create Command Figure
%  Main Figure: Initialize ++++++++++++++++++++++++++++++++++++++
% 'menubar','none' display tools: zoom
checksignals=figure('numbertitle','off',...
            'position',[46 42 600 450],...
            'keypressfcn',@manual_processing_ctrl,...
            'CloseRequestFcn',@close_and_update);
%% Video Switcher/SIGNAL NAVIGATOR
VideoSwitcher=uicontextmenu(checksignals);
checksignals.UIContextMenu=VideoSwitcher;
NextVideo=uimenu(VideoSwitcher,'Label','Go to','Callback',@gotosignal);
NextVideo=uimenu(VideoSwitcher,'Label','Next Video ->','Callback',@switchvideo);
PrevVideo=uimenu(VideoSwitcher,'Label','<- Previous Video','Callback',@switchvideo);
%% Define colormap for raster image: white/black:activity/orange:empty
CM=hot(20);         % Choose Colors
CM=CM(end:-1:1,:);  % Turn up-down
colormap(CM);       % Make Color map
%% Initialize AXIS and PLOTS (with randomw stuff) ##########################
ax1=subplot(3,2,1); % Raw Data - - - - - - - - - - - - - -      
plotsignal=plot(randn(1,10),'Parent',ax1);
plotsignal.Parent.YLabel.String='F_{raw}';
plotsignal.Parent.XTick=[];
hold(ax1,'on');
plotF0=plot(randn(1,10),'Parent',ax1);
plotF0.Color='m';
plotF0.LineStyle='-.';
plotF0.LineWidth=2;
% plotF0.Parent.XTick=[];
hold(ax1,'off');
ax2=subplot(3,2,3); % Detrended Data - - - - - - - - - - -
set(ax2,'NextPlot','replacechildren');
ax3=subplot(3,2,5); %  Driver - - - - - - - - - - - - - - - 
plotdriver=bar(randn(1,10),'Parent',ax3);
plotdriver.Parent.YLabel.String='Activity Driver';
plotdriver.Parent.XLabel.String='[min]';
plotdriver.Parent.XLabel.FontSize=9;
plotdriver.Parent.YLabel.FontSize=9;
plotdriver.Parent.TickLength=[0,0];
set(plotdriver,'ButtonDownFcn',@remove_driver)
ax4=subplot(3,2,[2,4]); % RASTER  - - - - - - - - - - - - - - - 
plotraster=imagesc(randn(10,10),'Parent',ax4);
plotraster.Parent.Title.String='RASTER';
plotraster.Parent.XTick=[];
plotraster.Parent.Box='off';
set(plotraster,'ButtonDownFcn',@remove_columns)
ax5=subplot(3,2,6); % CAG Signal  - - - - - - - - - - - - - - -
set(ax5,'NextPlot','replacechildren');
linkaxes([ax1,ax2,ax3],'x');
linkaxes([ax4,ax5],'x');
% Help Button
ButtonHelp=uicontrol('Parent',checksignals,'Style','pushbutton','String',...
    '?','Units','normalized','Position',[0.95 0.95 0.05 0.05],...
    'FontSize',10,'FontWeight','bold',...
    'Visible','on','Callback',@HelpCalciumMagic);

%% Initialize Plot **************************************************
n=1;
while isempty(indxSIGNALS{n})
    n=n+1;
    disp('... Searching Signal Indexes>>')
end
Nusefulpairs=size(IndxPair,1);
if n>Nusefulpairs
    n=1;
end
i=IndxPair(n,2); % CONDITIONS
j=IndxPair(n,1); % VIDEOS
[Cells,Frames]=size(SIGNALS{j,i});      % ALL NCells & Frames
isSIGNAL=indxSIGNALS{j,i};
% FR=Responses{j,i}(isSIGNAL,:);
% Here: Detect if signals are detected or Rejected      (!!!!)
if ~isempty(find(preLAMBDAS{j,i}(isSIGNAL)==0))
    RubricTitle=' Ca++ Signals | Condition: ';
    % Initialize Outputs
    for nvid=1:Npairs
        iaux=IndxPair(nvid,2); % Conditions
        jaux=IndxPair(nvid,1); % Videos
        indxSIGNALSOK{jaux,iaux}=setdiff(1:Cells,indxSIGNALS{jaux,iaux});
    end
    % undetectedindx=true;
else
    RubricTitle='Ca++ Fluorescence Signals | Condition: ';
    indxSIGNALSOK=indxSIGNALS;  % Detected INDEXES
    % undetectedindx=false;
end

checksignals.Name=[RubricTitle,' Exp ID: ',Experiment,' ',Names_Conditions{i},' Video: ',num2str(j)];
% Initialize Variables:
X=[]; XD=[]; D=[]; lambdass=[];
R=[]; X_SPARSE=[]; snrS=[];
FR=[]; 
% indxSIGNALSOK=[]; isSIGNALS=[];

retrieve_global_data();                 % Get Data from the Global Variables
indx_neuron=1;                          % Starting Index
neuron=isSIGNAL(indx_neuron);                   % Neuron
x=X(indx_neuron,:);                             % Raw Signal
r=FR(indx_neuron,:);                            % Response
xd=XD(indx_neuron,:);                           % Detrended
d=D(indx_neuron,:);                             % Driver
x_sparse=X_SPARSE(neuron,:);               % Clean Signal
snrc=snrS(neuron);                         % SNR
lambda=lambdass(indx_neuron);                   % Lambda of Preprocess
Get_Data();
Plot_Data_Now();
disp('... Ready!')
detectedneurons=[];
% isSIGNALS=[];
% ... and magic begins.
%% Nested Functions #######################################################
    function retrieve_global_data()
        % Read Data ++++++++++++++++++++++++++++++++++++++++++
        X=SIGNALS{j,i}(isSIGNAL,:);         % Raw Fluorescence {only-detected}
        XD=DETSIGNALS{j,i}(isSIGNAL,:);     % Detrended Fluorescence {only-detected}
        FR = Responses{j,i}(isSIGNAL,:);    % Response Funtions {only-detected}
        D=preDRIVE{j,i}(isSIGNAL,:);        % Driver Signals  {only-detected}
        R=RASTER{j,i};                      % Raster        {ALL CELLS}
        lambdass=preLAMBDAS{j,i}(isSIGNAL); % List of lambdas {only-detected}
        okReview=find(preLAMBDAS{j,i}(isSIGNAL)==0); % Check wich ones haven't been processed
        if ~isempty(okReview)
            disp('Searching in Undetected Ca++ Transients');
            hwait = figure('units','pixels',...
                'NumberTitle','off',...
                'position',[400 400 400 150],...
                'windowstyle','modal');
            
            uicontrol('style','text','string','>>Processing...',...
                'FontSize',24,...
                'BackgroundColor','black',...
                'ForegroundColor','green',...
                'units','pixels','position',[0 50 400 50]);
            drawnow;
            % Calculate Paramters:
            % Arbitrary Initiate to static workspace:
            p=1; L=1; taus_0=[1,1,1]; SpuriousTh=1; MaxLambdaMax=1;
            Load_Default_Values_SP;
            [FR(okReview,:),~,~]=AR_Estimation(XD(okReview,:),p,fs,L,taus_0);
            for k=1:length(isSIGNAL(okReview))
                if isempty(findpeaks( FR(okReview(k),:) ) )
                    FR(okReview(k),:)=-FR(okReview(k),:);
                    disp('WARNING: Response Function Missestimation')
                end
            end
            Responses{j,i}(isSIGNAL,:)=FR;
            % Input:        empty-> +&- | 0->only+ | 1->+or-
            [D(okReview,:),lambdass(okReview)]=maxlambda_finder(XD(okReview,:),FR(okReview,:),1);  % [[[DECONVOLUTION]]]
            % Get Raster: Driver Mode
            for k=1:numel(okReview)
                switch RasterAlgorithm
                case 'Driver'
                    R(isSIGNAL(okReview(k)),:)=get_raster(1,D(okReview(k),:),1);
                case 'Derivative'
                    R(isSIGNAL(okReview(k)),:)=get_raster(3,D(okReview(k),:),1,FR(okReview(k),:));
                end
                % R(isSIGNAL(okReview(k)),D(okReview(k),:)>0)=1;
            end
            preLAMBDAS{j,i}(isSIGNAL(okReview))=lambdass(okReview);
            preDRIVE{j,i}(isSIGNAL(okReview),:)=D(okReview,:);
            close(hwait);
            disp('Done.');
            % Plot_Data_Now;
        end
        
        if isempty(SIGNALSclean{j,i})
            fprintf('>> Getting Sparse Signal...')
            SIGNALSclean{j,i}=zeros(Cells,Frames);
            SNRlambda{j,i}=zeros(Cells,1);
            X_SPARSE=zeros(Cells,Frames);        % Initialize Sparse Clean Signals
            snrS=zeros(Cells,1);     % Initialize Output SNRs
            for k=1:length(isSIGNAL)        % {only-detected}
                xd=XD(k,:);
                d=D(k,:);
                r=FR(k,:);
                x_sparse=sparse_convolution(d,r);
                X_SPARSE(isSIGNAL(k),:)=x_sparse;     % Sparse Signals
                snrc=10*log(var(x_sparse)/var(xd'-x_sparse));
                snrS(isSIGNAL(k))=snrc;               % SNR Sparse
                fprintf('.')
            end
            fprintf(' Done\n')
        else
            fprintf('>> Getting Sparse Signal...')
            X_SPARSE=SIGNALSclean{j,i};
            % PATCH SNRlambda shoul've been saved
            if ~isempty(SNRlambda{j,i})
                snrS=SNRlambda{j,i};
            else
                snrS=10*log(var(XD')./var(X_SPARSE(isSIGNAL,:)'));
                SNRlambda{j,i}=snrS;
                disp('~fix~');
            end
            fprintf(' Already Done\n')
        end
    end

    function Get_Data()
        neuron=isSIGNAL(indx_neuron);         % Neuron
        x=X(indx_neuron,:);                   % Raw Signal
        r=FR(indx_neuron,:);                  % Response
        xd=XD(indx_neuron,:);                 % Detrended
        d=D(indx_neuron,:);                   % Driver
        d(d<0)=0;
        x_sparse=X_SPARSE(neuron,:);     % Clean Signal
        snrc=10*log(var(x_sparse)/var(xd-x_sparse));   % SNR
        lambda=lambdass(indx_neuron);  
    end

    function Plot_Data_Now()
        % Raw Data - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        plotsignal.YData=x;
        plotF0.YData=x-xd;
        plotsignal.Parent.XLim=[1,Frames];
        plotsignal.Parent.YLim=[min(x),max(x)];
        plotsignal.Parent.XGrid='on';
        plotsignal.Parent.YGrid='on';
        title(plotsignal.Parent,['Neuron: ',num2str(neuron),...
            ' Progress: ',num2str(100*indx_neuron/length(isSIGNAL),3),'%',],...
              'FontSize',8,'FontWeight','normal');
        % Detrended Data - - - - - - - - - - - - - - - - - - - - - - - - - 
        cla(ax2);
        hold(ax2,'on')
        plot(ax2,xd,'b','LineWidth',1)
        grid(ax2,'on');
        plot(ax2,x_sparse,'g','LineWidth',2)
        plot(ax2,[0,numel(xd)],[std(xd-x_sparse),std(xd-x_sparse)],'-.r');
        plot(ax2,[0,numel(xd)],-[std(xd-x_sparse),std(xd-x_sparse)],'-.r');
        hold(ax2,'off')
        axis(ax2,'tight');
        ax2.XTick=0:30*fs:ceil(Frames/(60*fs))*60*fs;
        ax2.XTickLabel=[];
        ylabel(ax2,'$\Delta F / F_0$','Interpreter','Latex')
        title(ax2,['SNR=',num2str(snrc),'[dB]'],...
              'FontSize',8,'FontWeight','normal');
        %  Driver - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        if max(d)<=0
            dplot=d;
            plotdriver.EdgeColor='r';
            plotdriver.LineWidth=2;
        else
            dplot=d/max(d);    % Normalize
            plotdriver.EdgeColor='k';
            plotdriver.LineWidth=0.5;
        end
        plotdriver.YData=dplot;
        axis(plotdriver.Parent,'tight');
        grid(plotdriver.Parent,'on');
        title(plotdriver.Parent,['\lambda=',num2str(lambda)],...
             'FontSize',8,'FontWeight','normal');
        plotdriver.Parent.XTick=0:60*fs:ceil(Frames/(60*fs))*60*fs;
        % Label Xticks:
        for ilab=1:numel(plotdriver.Parent.XTick)
            plotdriver.Parent.XTickLabel{ilab}=num2str(ilab-1);
        end
        % RASTER  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Raux=R*2;
        if ~isempty(Raux(neuron,:)>0)
            Raux(neuron,Raux(neuron,:)>0)=1;
        end
        plotraster.CData=Raux;
        plotraster.Parent.YLim=[0,Cells+1];
        plotraster.Parent.XLim=[1,Frames];
        plotraster.Parent.XTick=[];
        plotraster.Parent.YTick=neuron;
        % CoAc Signal  - - - - - - - - - - - - - - - - - - - - - - - - - - 
        CoAcGraphy=sum(R);
        cla(ax5); % Replot Stuff
        hold(ax5,'on')
        plot(ax5,CoAcGraphy,'b','LineWidth',2)
        plot(ax5,find(d>0),CoAcGraphy(d>0),'LineStyle','none','Marker','.','MarkerSize',10,'Color',[0.6,0,0.6])
        grid(ax5,'on');
        hold(ax5,'off')
        axis(ax5,'tight');
        ylabel(ax5,'CAG')
        ax5.XTick=0:30*fs:ceil(Frames/(60*fs))*60*fs;
        ax5.XTickLabel=[];
        drawnow;
    end

    function manual_processing_ctrl(checksignals,~,~)
        key=get(checksignals,'CurrentKey');
        if strcmp(key,'rightarrow')
            disp('Next ------->')
            if indx_neuron<length(isSIGNAL)
                indx_neuron=indx_neuron+1; % Increases
            else
                indx_neuron=length(isSIGNAL); % Stays
            end
            Get_Data;
            disp(neuron)
            update_data;
            Plot_Data_Now;
            % SHOULD SAVE DATA
        elseif strcmp(key,'leftarrow')
            disp('<-- Previous')
            if indx_neuron==1
                indx_neuron=1;
            else
                indx_neuron=indx_neuron-1;
            end
            Get_Data;
            disp(neuron)
            update_data;
            Plot_Data_Now;
            % Set_Name_Axis
        elseif strcmp(key,'uparrow')
            disp('Increase lambda /\')
            disp('............... |')
            Get_Data;
            lambda=(1+delta_lambda)*lambda;
            reprocess_data; 
            Plot_Data_Now;
            update_data;
        elseif strcmp(key,'downarrow')
            disp('...............  |')
            disp('Decrease lambda \/')
            Get_Data;
            lambda=(1-delta_lambda)*lambda;
            reprocess_data;
            Plot_Data_Now;
            update_data;
        end
    end


    function gotosignal(~,~)
        GoToN=inputdlg('ROI number');
        if ~isempty(GoToN)
            Ngosignal=round(str2double(GoToN));
            if sum(ismember(isSIGNAL,Ngosignal))>0
                fprintf('>Moving to %i ROI\n',Ngosignal);
                indx_neuron=find(isSIGNAL==Ngosignal);
            else
                [~,indx_neuron]=min(abs(isSIGNAL-Ngosignal));
                fprintf('>No signal for %i ROI, but %i is the closest\n',Ngosignal,isSIGNAL(indx_neuron));
            end
            Get_Data;
            disp(neuron)
            update_data;
            Plot_Data_Now;
        else
            fprintf('>No ROI selected\n');
        end
    end

    function reprocess_data()
        [d,x_sparse,~]=magic_sparse_deconvolution(xd,r,lambda);
        % re-re-process if detrending issue: emtpy last input
        % [d,~,~,~,~,~]=analyze_driver_signal(d',r,xd,x_sparse);
        % Onky anlyzes Driver: any last input
        if ~isempty(d(d>0))
            [d,~,x_sparse,~,~,~]=analyze_driver_signal(d',r,xd,x_sparse,0);
        end
        switch RasterAlgorithm
            case 'Driver'
                R(neuron,:)=get_raster(1,d,1);
            case 'Derivative'
                R(neuron,:)=get_raster(3,d,1,FR(indx_neuron,:));
        end
        % By Driver
        % R(neuron,:)=0;
        %R(neuron,d>0)=1;
        % By (+) Derivative:
        % R(neuron,diff(x_sparse)>0)=1;
        if size(xd,1)~=size(x_sparse,1)
            x_sparse=x_sparse';
        end
        sparsenoise=xd-x_sparse;
        snrc=10*log(var(x_sparse)/var(sparsenoise)); 
    end

    function update_data()
        % XD(indx_neuron,:)=xd;          
        D(indx_neuron,:)=d;
        snrS(neuron)=snrc;
        X_SPARSE(neuron,:)=x_sparse;
        lambdass(indx_neuron)=lambda;
        switch RasterAlgorithm
            case 'Driver'
                R(neuron,:)=get_raster(1,d,1);
            case 'Derivative'
                R(neuron,:)=get_raster(3,d,1,FR(indx_neuron,:));
        end
        % Update Raster: Driver Mode
        % R(neuron,:)=0;
        % R(neuron,d>0)=1;
        % Derivative Mode
        % R(neuron,diff(x_sparse)>0)=1;
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
        d(frames2del(1):frames2del(2))=0;
        % x_sparse(frames2del(1):frames2del(2))=0;
        update_data;        % saves daux
        Plot_Data_Now;      % update plot
        figure(checksignals) % Focus on Window
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
        D(:,frames2del(1):frames2del(2))=0;
        Plot_Data_Now;      % update plot
        % update_data;        % saves daux
        figure(checksignals) % Focus on Window
        sprintf('Deleted Raster %d : %d frames \n',frames2del)
    end 
 
    function switchvideo(source,~)
        fprintf('Getting Data...\n')
        update_global_data;
        switch source.Label
            case 'Next Video ->'                
                if n+1>Npairs
                    n=Npairs;
                    i=IndxPair(n,2); % CONDITIONS
                    j=IndxPair(n,1); % VIDEOS
                    msgbox({'END';'OF';'THE';'ROAD'});
                    disp('END')
                else
                    n=n+1;
                    i=IndxPair(n,2); % CONDITIONS
                    j=IndxPair(n,1); % VIDEOS
                    disp('               Next Video=>')
                end
            case '<- Previous Video'
                if n-1<1
                    n=1;
                    i=IndxPair(n,2); % CONDITIONS
                    j=IndxPair(n,1); % VIDEOS
                    msgbox({'BEGINNING';'OF';'THE';'ROAD'});
                    disp('END')
                else
                    n=n-1;
                    i=IndxPair(n,2); % CONDITIONS
                    j=IndxPair(n,1); % VIDEOS
                    disp('                  <=Previus Video')
                end
        end
        checksignals.Name=[RubricTitle,'Exp ID: ',Experiment,' ',Names_Conditions{i},' Video: ',num2str(j)];
        indx_neuron=1; 
        isSIGNAL=indxSIGNALS{j,i};
        [Cells,Frames]=size(SIGNALS{j,i});
        retrieve_global_data;
        Get_Data;
        Plot_Data_Now;
        fprintf('... Data Ready\n')
    end

    function update_global_data()
        % Read Data ++++++++++++++++++++++++++++++++++++++++++
        detectedneurons=find(sum(R,2)>0);       % Index of Actual Detected Neurons        
        indxSIGNALSOK{j,i}=detectedneurons;     % Update ACtual Detected Ca++Trans
        preDRIVE{j,i}(isSIGNAL,:)=D;            % Driver Signals  {only-detected}
        preLAMBDAS{j,i}(isSIGNAL)=lambdass;     % List of lambdas {only-detected}
        RASTER{j,i}=R;                          % Raster        {ALL CELLS}
        SIGNALSclean{j,i}=X_SPARSE;             % Clean Sparse Signals {ALL CELLS}
        SNRlambda{j,i}=snrS;          % Signal Noise Ratio {only-detected}
        %if undetectedindx
        %    Responses{j,i}(isSIGNAL)=FR;
        %end
        disp('Data: [Updated]')
    end

    function close_and_update(~,~,~)
        update_global_data;
        uiwait(msgbox('Open .mat file to overwrite & update changes'));
        % Save as NEW Cell of Checked Indexes 
        % isSIGNALS=indxSIGNALSOut; % To Save
        % Update NO DETECTDED SIGNALS
        % for nvid2=1:Npairs
        %    i=IndxPair(nvid2,2); % Conditions
        %    j=IndxPair(nvid2,1); % Videos
        %    notSIGNALS{j,i}=setdiff(1:Cells,isSIGNALS{j,i});
        % end
        checkname=1;
        while checkname==1
            % Get Directory
            DP=pwd;
            Slashes=find(DP==filesep);
            fprintf('\n\n [ To update changes, please select the .mat file ] \n\n')
            DefaultPath=fullfile(DP(1:Slashes(end)),DirectorySelector(1));
            if ~isdir(DefaultPath)
                DefaultPath=pwd; % Current Diretory of FinderSpiker
            end
            [FileName,PathName] = uigetfile('*.mat',[' Pick the Analysis File to SAVE CHANGES ',Experiment],...
                'MultiSelect', 'off',DefaultPath);
            % dotindex=find(FileName=='.');
            if strcmp(FileName(1:end-4),Experiment)
                checkname=0;
                % SAVE DATA
                save([PathName,FileName],'preDRIVE','preLAMBDAS','RASTER','RasterAlgorithm',...
                    'Responses','SIGNALSclean','SNRlambda','indxSIGNALSOK','-append');
                disp([Experiment,'   -> RASTERS UPDATED (Visual ~ Inspection)'])
                delete(gcf)
            elseif FileName==0
                checkname=0;
                disp('....CANCELLED')
                delete(gcf)
            else
                disp('Not the same Experiment!')
                disp('Try again!')
            end
        end    
        delete(checksignals);
    end
    function HelpCalciumMagic(~,~,~)
        HelpGUI;
    end
end
%% END OF THE WORLD