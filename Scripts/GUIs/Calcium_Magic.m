% Function to Clean Manuallly Calcium Transients
% Input
% indxSIGNAL: Cell of Indexed Signal as Detected Ca^2+ Transient
% Global Inputs:
%   SIGNALS             Raw Fluuorescence
%   DETSIGNALS          Detrended FLuorescence
%   preDRIVE            Preliminar Drive
%   preLAMBDAS          Preliminar Lambdas
%   RASTER              Detected Raster ACtivity
%   Responses           Estimated Dye Response
%   Names_Conditions    Names of Conditions
% Outputs
%   indxSINGALSOK:      Acutal Detected Signals
%   SIGNALSclean:       Sparse Cleaned up Signals
%   SNRs                Signal Noise Ratio
function indxSINGALSOK = Calcium_Magic(indxSIGNALS)
%% Setup
disp('Loading...')
% Call Global Variables
global SIGNALS;             % Raw Fluorescence Signals
global DETSIGNALS;          % Detrended Fluorescence Signals
global preDRIVE;            % Automatic Founded Driver Function
global preLAMBDAS;          % Automatic Lambdas Founded
global RASTER;              % Raster Activity
global Responses;           % Fluorescence Response Signals
global Names_Conditions;    % Names of The Conditions
global SIGNALSclean;        % Sparse Clean Signals
global SNRlambda;           % Signal Noise Ratio by Sparse Deconvolution

% Initialize Outputs
indxSINGALSOK=indxSIGNALS;  % Detected INDEXES
% Setup Variables
delta_lambda=0.25;          % \lambda Step lambda_next=(1+_deltalambda)*lambda_present
[NVmax,NC]=size(SIGNALS);   % NC: N Conditions
% Get Useful Pair of Indexes i,j
IndxPair=[];
for i=1:NC
    for j=1:NVmax
       if ~isempty(SIGNALS{j,i})
           IndxPair=[IndxPair;j,i];
       end
    end
end
Npairs=size(IndxPair,1);

%% Create Command Figure
%  Main Figure: Initialize ++++++++++++++++++++++++++++++++++++++
% 'menubar','none' display tools: zoom
checksignals=figure('numbertitle','off',...
            'position',[46 42 600 450],...
            'keypressfcn',@manual_processing_ctrl);
% Video Switcher
VideoSwitcher=uicontextmenu(checksignals);
checksignals.UIContextMenu=VideoSwitcher;
NextVideo=uimenu(VideoSwitcher,'Label','Next Video ->','Callback',@switchvideo);
PrevVideo=uimenu(VideoSwitcher,'Label','<- Previous Video','Callback',@switchvideo);
% Define colormap for raster image: white/black:activity/orange:empty
CM=hot(20);         % Choose Colors
CM=CM(end:-1:1,:);  % Turn up-down
colormap(CM);       % Make Color map
% Initialize AXIS and PLOTS (with randomw stuff) ##########################
ax1=subplot(3,2,1); % Raw Data - - - - - - - - - - - - - -      
plotsignal=plot(randn(1,10),'Parent',ax1);
plotsignal.Parent.YLabel.String='F_{raw}';
plotsignal.Parent.XTick=[];
ax2=subplot(3,2,3); % Detrended Data - - - - - - - - - - -
set(ax2,'NextPlot','replacechildren');
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
ax5=subplot(3,2,6); % CAG Signal  - - - - - - - - - - - - - - -
set(ax5,'NextPlot','replacechildren');
linkaxes([ax1,ax2,ax3],'x');
linkaxes([ax4,ax5],'x');
%% Initialize Plot **************************************************
n=1;
i=IndxPair(n,2); % CONDITIONS
j=IndxPair(n,1); % VIDEOS

isSIGNAL=indxSIGNALS{j,i};
% Here: Detect if signals are detected or Rejected      (!!!!)
RubricTitle='Detected Ca++ Transients | Condition: ';
[Cells,Frames]=size(SIGNALS{j,i});      % ALL NCells & Frames
checksignals.Name=[RubricTitle,Names_Conditions{i},' Vid: ',num2str(j)];
% Initialize Variables:
X=[]; FR=[]; XD=[]; D=[]; lambdass=[];
R=[]; X_SPARSE=[]; snrS=[];

retrieve_global_data();
indx_neuron=1; 
neuron=isSIGNAL(indx_neuron);                   % Neuron
x=X(indx_neuron,:);                             % Raw Signal
r=FR(indx_neuron,:);                            % Response
xd=XD(indx_neuron,:);                           % Detrended
d=D(indx_neuron,:);                             % Driver
x_sparse=X_SPARSE(indx_neuron,:);               % Clean Signal
snrc=snrS(indx_neuron);                         % SNR
lambda=lambdass(indx_neuron);                   % Lambda of Preprocess
Get_Data();
Plot_Data_Now();
disp('... Ready!')
% ... and magic begins.
%% Nested Functions #######################################################
    function retrieve_global_data()
        % Read Data ++++++++++++++++++++++++++++++++++++++++++
        X=SIGNALS{j,i}(isSIGNAL,:);         % Raw Fluorescence {only-detected}
        FR = Responses{j,i}(isSIGNAL,:);    % Response Funtions {only-detected}
        XD=DETSIGNALS{j,i}(isSIGNAL,:);     % Detrended Fluorescence {only-detected}
        D=preDRIVE{j,i}(isSIGNAL,:);        % Driver Signals  {only-detected}
        lambdass=preLAMBDAS{j,i}(isSIGNAL); % List of lambdas {only-detected}
        R=RASTER{j,i};                      % Raster        {ALL CELLS}
        if isempty(SIGNALSclean{j,i})
            fprintf('>> Getting Sparse Signal...')
            SIGNALSclean{j,i}=zeros(size(SIGNALS{j,i}));
            SNRlambda{j,i}=zeros(size(SIGNALS{j,i},1),1);
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
            fprintf(' Done\n')
        else
            fprintf('>> Getting Sparse Signal...')
            X_SPARSE=SIGNALSclean{j,i}(isSIGNAL,:);
            snrS=SNRlambda{j,i}(isSIGNAL);
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
        x_sparse=X_SPARSE(indx_neuron,:);     % Clean Signal
        snrc=10*log(var(x_sparse)/var(xd-x_sparse));   % SNR
        lambda=lambdass(indx_neuron);  
    end

    function Plot_Data_Now()
        % Raw Data - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        plotsignal.YData=x;
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
        ylabel(ax2,'$\Delta F / F_0$','Interpreter','Latex')
        title(ax2,['SNR=',num2str(snrc),'[dB]'],...
              'FontSize',8,'FontWeight','normal');
        %  Driver - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        dplot=d/max(d);    % Normalize
        plotdriver.YData=dplot;
        axis(plotdriver.Parent,'tight');
        grid(plotdriver.Parent,'on');
        title(plotdriver.Parent,['\lambda=',num2str(lambda)],...
             'FontSize',8,'FontWeight','normal');
        % RASTER  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Raux=R*2;
        if ~isempty(Raux(neuron,:)>0)
            Raux(neuron,Raux(neuron,:)>0)=1;
        end
        plotraster.CData=Raux;
        plotraster.Parent.YLim=[0,Cells];
        plotraster.Parent.XLim=[1,Frames];
        plotraster.Parent.XTick=[];
        plotraster.Parent.YTick=neuron;
        % CoAc Signal  - - - - - - - - - - - - - - - - - - - - - - - - - - 
        CoAcGraphy=sum(R);
        cla(ax5)
        hold(ax5,'on')
        plot(ax5,CoAcGraphy,'b','LineWidth',2)
        plot(ax5,find(d>0),CoAcGraphy(d>0),'LineStyle','none','Marker','.','MarkerSize',10,'Color',[0.6,0,0.6])
        grid(ax5,'on');
        hold(ax5,'off')
        axis(ax5,'tight');
        ylabel(ax5,'CAG')
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
    
    function reprocess_data()
        [d,x_sparse,~]=magic_sparse_deconvolution(xd,r,lambda);
        [d,~,~,~,~,~]=analyze_driver_signal(d',r,xd,x_sparse);
        % By Driver
        R(neuron,:)=0;
        R(neuron,d>0)=1;
        % By (+) Derivative:
        % R(neuron,diff(x_sparse)>0)=1;
        x_sparse=x_sparse';
        sparsenoise=xd-x_sparse;
        snrc=10*log(var(x_sparse)/var(sparsenoise)); 
    end

    function update_data()
        XD(indx_neuron,:)=xd;          
        D(indx_neuron,:)=d;
        snrS(indx_neuron)=snrc;
        X_SPARSE(indx_neuron,:)=x_sparse;
        lambdass(indx_neuron)=lambda;
        % Update Raster: Driver Mode
        R(neuron,:)=0;
        R(neuron,d>0)=1;
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
        checksignals.Name=[RubricTitle,Names_Conditions{i},' Vid: ',num2str(j)];
        indx_neuron=1; 
        isSIGNAL=indxSIGNALS{j,i};
        retrieve_global_data;
        Get_Data;
        Plot_Data_Now;
        fprintf('... Data Ready\n')
    end

    function update_global_data()
        % Read Data ++++++++++++++++++++++++++++++++++++++++++
        detectedneurons=find(sum(R,2)>0);       % Index of Actual Detected Neurons        
        indxSINGALSOK{j,i}=detectedneurons;     % Update ACtual Detected Ca++Trans
        preDRIVE{j,i}(isSIGNAL,:)=D;            % Driver Signals  {only-detected}
        preLAMBDAS{j,i}(isSIGNAL)=lambdass;     % List of lambdas {only-detected}
        RASTER{j,i}=R;                          % Raster        {ALL CELLS}
        SIGNALSclean{j,i}(isSIGNAL,:)=X_SPARSE; % Clean Sparse Signals {ALL CELLS}
        SNRlambda{j,i}(isSIGNAL)=snrS;          % Signal Noise Ratio {only-detected}
        disp('Data: [Updated]')
    end
%% Exit to Update
end
%% END OF THE WORLD