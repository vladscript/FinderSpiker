function varargout = NeuralNetwork(varargin)
% NEURALNETWORK M-file for NeuralNetwork.fig
%      NEURALNETWORK, by itself, creates a new NEURALNETWORK or raises the existing
%      singleton*.
%
%      H = NEURALNETWORK returns the handle to a new NEURALNETWORK or the handle to
%      the existing singleton*.
%
%      NEURALNETWORK('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEURALNETWORK.M with the given input arguments.
%
%      NEURALNETWORK('Property','Value',...) creates a new NEURALNETWORK or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NeuralNetwork_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NeuralNetwork_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help NeuralNetwork
% Last Modified by GUIDE v2.5 18-Sep-2014 11:09:16

%% ...::: Initialization :::...
    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @NeuralNetwork_OpeningFcn, ...
                       'gui_OutputFcn',  @NeuralNetwork_OutputFcn, ...
                       'gui_LayoutFcn',  [] , ...
                       'gui_Callback',   []);
    if nargin && ischar(varargin{1})
        gui_State.gui_Callback = str2func(varargin{1});
    end
    if nargout
        [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
    else
        gui_mainfcn(gui_State, varargin{:});
    end
    % End initialization code - DO NOT EDIT
end
% --- Executes just before NeuralNetwork is made visible.
function NeuralNetwork_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject;
    global gray_image yellow_image green_image waiting_image
    global step_execution 
    gray_image = imread('execute.png');
    yellow_image = imread('warning.png');
    green_image = imread('done.png');
    waiting_image = imread('waiting.png');
    set(handles.PlotDataPushbutton,'CData', gray_image);
    set(handles.PlotPeaksPushbutton,'CData', gray_image);
    set(handles.DimRedPushbutton,'CData', gray_image);
    set(handles.ClusteringPushbutton,'CData', gray_image);
    set(handles.PlotNetworkPushbutton,'CData', gray_image);
    
    step_execution = [0 0 0 0 0];
    guidata(hObject, handles);
end
function varargout = NeuralNetwork_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;
end

%% ...::: Functions :::...
%% Processing
function Processing(mode,message,handles)
    switch mode
        case 'msj'
            set(handles.ProcessText,'string',strcat('Processing... ',message))
        case 'error'
            set(handles.ProcessText,'string',strcat('Error!... ',message))
        case 'end'
            set(handles.ProcessText,'string',strcat('Done!... ',message))
    end
end
%% Get Experiment data
function Experiment=Get_Experiment(handles)
    DataNum=get(handles.DataPopupmenu,'Value');
    DataStrings=get(handles.DataPopupmenu,'String');
    DataName=DataStrings{DataNum};
    if evalin('base',['exist(''' DataName '_Analysis'')']);
        Experiment=evalin('base',[DataName '_Analysis']);
    else
        Experiment=[];
    end
end

%% Get data
function Experiment=Get_Data(handles)
    % Get data
    DataNum=get(handles.DataPopupmenu,'Value');
    if DataNum==1
        Data=[];
        DataName='';
    else
        DataStrings=get(handles.DataPopupmenu,'String');
        DataName=DataStrings{DataNum};
        Data=evalin('base',['exist(''' DataName ''')']);
        if Data
            Data=evalin('base',DataName);
        else
            Data=[];
            DataName='';
            UpdateWorkspaceMenu_Callback([], [], handles);
        end
            
    end
    % Evaluate data
    %set(handles.InfoText,'ForegroundColor','black')
    [F C]=size(Data);
    if (length(size(Data))==2 && F>1 && C>1)
        fps=str2num(get(handles.FPSEdit,'string'));
        t=F/fps;
        if (F<500||C>300)
            Warn='Warning: may not be an experiment.\n\n';
            %set(handles.InfoText,'ForegroundColor','red')
        else
            Warn='';
        end
        info=sprintf([Warn DataName '\n\n %d neurons\n %d frames\n'...
            '%0.1f seconds'],C,F,t);
    else
        info='Warning: select a valid experiment.';
        %set(handles.InfoText,'ForegroundColor','red')
        Data=[];
        DataName='';
    end
    %set(handles.InfoText,'String',info)
    
    [Frames Cells]=size(Data);
    
    % Outputs
    Experiment.Data.Data=Data;
    Experiment.Data.DataName=DataName;
    Experiment.Data.Neurons=Cells;
    Experiment.Data.Frames=Frames;
    
end

%% Get XY
function Experiment=Get_XY(Experiment,handles)
    % Get data
    XYNum=get(handles.XYPopupmenu,'Value');
    if XYNum==1
        XY=[];
    else
        XYStrings=get(handles.XYPopupmenu,'String');
        XYName=XYStrings{XYNum};
        XY=evalin('base',['exist(''' XYName ''')']);
        if XY
            XY=evalin('base',XYName);
        else
            XY=[];
            UpdateWorkspaceMenu_Callback([], [], handles);
        end
    end
    
    % Evaluate data
    set(handles.XYPopupmenu,'BackgroundColor','white')
    C=Experiment.Data.Neurons;
    if ~((length(size(XY))==2 && size(XY,2)==2 && size(XY,1)==C) || isempty(XY))
        set(handles.XYPopupmenu,'BackgroundColor','red')
        XY=[];
    end
    Experiment.Data.XY=XY;
end

%% Get XY colors
function Experiment=Get_XYColors(Experiment,handles)

    C=Experiment.Data.Neurons;
    XYNum=get(handles.XYColorsPopupmenu,'Value');
    set(handles.XYColorsPopupmenu,'BackgroundColor','White')
    if XYNum==1
        XYColors=repmat([0 0 0],C,1);
    else
        XYStrings=get(handles.XYColorsPopupmenu,'String');
        XYName=XYStrings{XYNum};
        XYColors=evalin('base',['exist(''' XYName ''')']);
        if XYColors
            XYColors=evalin('base',XYName);
        else
            XYColors=[];
            UpdateWorkspaceMenu_Callback([], [], handles);
        end

        [NCellColor NColors]=size(XYColors);
        if (NColors==3 && max(max(XYColors))<=1 && min(min(XYColors))>=0)
            if NCellColor==1
                XYColors=repmat(XYColors,C,1);
            elseif NCellColor~=C
                XYColors=repmat([0 0 0],C,1); % color: black
                set(handles.XYColorsPopupmenu,'BackgroundColor','Red')
                %errordlg('Dimensions are not consistent.','Error')
            end
        else
            XYColors=repmat([0 0 0],C,1); % color: black
            set(handles.XYColorsPopupmenu,'BackgroundColor','Red')
            %errordlg('Dimensions are not consistent.','Error')
        end
    end
    
    % Evaluate XY colors
    %{
    [NCellColor NColors]=size(XYColors);
    if (NColors==3 && max(max(XYColors))<=1 && min(min(XYColors))>=0)
        if NCellColor==1
            XYColors=repmat(XYColors,C,1);
        elseif NCellColor~=C
            XYColors=repmat([0 0 0],C,1); % color: black
        end
    else
        XYColors=repmat([0 0 0],C,1); % color: black
    end
    %}
    
    
    
    Experiment.Colors.XY=XYColors;
end

%% Get group colors
function Experiment=Get_GroupColors(Experiment,handles)
    GroupNum=get(handles.GroupColorsPopupmenu,'Value');
    set(handles.GroupColorsPopupmenu,'BackgroundColor','White')
    default=[0.9 0.0 0.0;...     % red
             0.5 0.8 0.0;...     % green
             0.0 0.5 1.0;...     % sky blue
             0.9 0.9 0.0;...     % yellow
             0.9 0.0 0.9;...     % magenta
             0.0 0.9 0.9;...     % cyan
             1.0 0.5 0.0;...     % orange 
             0.5 0.0 1.0;...     % purple
             1.0 0.7 0.7;...     % rose
             0.6 0.3 0.0];       % brown
    if GroupNum==1
        GroupColors=default;
    else
        GroupStrings=get(handles.GroupColorsPopupmenu,'String');
        GroupName=GroupStrings{GroupNum};
        GroupColors=evalin('base',['exist(''' GroupName ''')']);
        if GroupColors
            GroupColors=evalin('base',GroupName);
            NColors=size(GroupColors,2);
            if (NColors==3 && max(max(GroupColors))<=1 && min(min(GroupColors))>=0)
                GroupColors=repmat(GroupColors,10,1);
            else
                GroupColors=default;
                set(handles.GroupColorsPopupmenu,'BackgroundColor','Red')
            end
        else
            GroupColors=default;
            UpdateWorkspaceMenu_Callback([], [], handles);
        end
        
    end
    Experiment.Colors.States=GroupColors;
end

%% Get States of Neurons
function Experiment=GetStatesNeuron(Experiment)

    % Inputs
    Data=Experiment.Peaks.ActiveDataPeaks;
    Idx=Experiment.Clustering.VectorStateIndex;
    Groups=Experiment.Clustering.TotalStates;
    
    
    Cells=size(Data,2);
    CellsGroups=zeros(Groups,Cells);
    CellsAcGroups=zeros(Groups,Cells);
    for i=1:Groups
        PeaksGroupIdx= Idx==i;
        CellsAcGroups(i,:)=sum(Data(PeaksGroupIdx,:),1);
        CellsGroup= logical(CellsAcGroups(i,:));
        CellsGroups(i,CellsGroup)=i;
    end
    
    % Outputs
    Experiment.Structure.NeuronsStates=CellsGroups;
    Experiment.Structure.RateNeuronsStates=CellsAcGroups;
    
end

%% Get States Combinations of Neuron
function Experiment=GetStatesComb(Experiment)
    
    % Inputs
    CellsGroups=Experiment.Structure.NeuronsStates;
    
    [Groups C]=size(CellsGroups);
    States=[]; Labels={}; Idx=[]; Indexes=[];
    Combinations=2^Groups-1;
    for i=1:Combinations
        bin=dec2bin(i,Groups);
        GroupsToCompare=[];
        Label=[];
        for j=1:Groups
            Compare=str2num(bin(j));
            if Compare
                GroupsToCompare=[GroupsToCompare; CellsGroups(j,:)>=1];
                if Label
                    Label=[Label ',' num2str(j)];
                else
                    Label=num2str(j);
                end
            else
                GroupsToCompare=[GroupsToCompare; CellsGroups(j,:)==0];
            end
        end
        Idx{i}=find(sum(GroupsToCompare,1)==Groups);
        States(i)=length(Idx{i})/C;
        Labels(i)={Label};
    end
    [Labels idxSort]=sort(Labels);
    States=States(idxSort);
    for i=1:Combinations
        Indexes=[Indexes Idx{idxSort(i)}];
        IndexesComb(Idx{i})=i;
    end
    
    % Outputs
    Experiment.Structure.StatesRatio=States;
    Experiment.Structure.LabelsStates=Labels;
    Experiment.Structure.Indexes=Indexes;
    Experiment.Structure.StatesCombinations=IndexesComb;
    Experiment.Structure.IndexesLabelsSort=idxSort;
    
end

%% Get adjacency matrix from peaks
function AdjacentPeaks=GetAdjacencyFromPeaks(handles,Data)
    [peaks C]=size(Data);
    AdjacentPeaks=zeros(C);
    
    CorrMethods=get(handles.CorrMethodPopupmenu,'string');
    idx=get(handles.CorrMethodPopupmenu,'value');
    CorrMethod=CorrMethods{idx};
    if strcmp(CorrMethod,'Synchronization')
        for i=1:peaks
            for j=1:C-1
                for k=j+1:C
                    x=Data(i,j);
                    y=Data(i,k);
                    if x&&y
                        AdjacentPeaks(j,k)=AdjacentPeaks(j,k)+1;
                        AdjacentPeaks(k,j)=AdjacentPeaks(k,j)+1;
                    end
                end
            end
        end
        %AdjacentPeaks=(AdjacentPeaks/peaks);
    elseif strcmp(CorrMethod,'Jaccard')
        AdjacentPeaks=squareform(1-pdist(Data','jaccard'));
    else
        AdjacentPeaks=corr(Data,'type',CorrMethod);
    end
    
end    

%% Get maximum coactivity and maximum activity
function [AcMax CoMax]=GetAcCoMax(Data)
    F=size(Data,1);
    Co=sum(Data,2);
    CoMax=max(Co);
    Ac=sum(Data,1)*100/F;
    AcMax=round(max(Ac));
end

%% Join Peaks
function Peaks = Join_Peaks(StatesIdx,PeaksIdx)
    
    NS=length(StatesIdx);
    NP=length(PeaksIdx);
    if NS~=NP
        for i=1:NS
            idx(i)=find(PeaksIdx==i,1,'first');
        end
    else
        idx=find(PeaksIdx);
    end
    
    F=length(idx);
    Peaks(1)=StatesIdx(1);
    j=2;
    for i=2:F
        if idx(i-1)~=idx(i)-1
            Peaks(j)=StatesIdx(i);
            j=j+1;
        end
    end
end

%% States transitions
function Transitions = StatesTrans(StatesSequence)
    peaks=length(StatesSequence);
    states=max(StatesSequence);
    Transitions=zeros(states);
    for i=1:states
        for j=1:states
            for k=1:peaks-1
                if StatesSequence(k)==i && StatesSequence(k+1)==j
                    Transitions(i,j)=Transitions(i,j)+1;
                end
            end
        end
    end
    Transitions=Transitions/(peaks-1);
end

%% Plot Arrow
function [x y] = PlotArrow(XY,Radius,LengthArrow,LineWidth,Color)

    % Curve
    l=norm(XY(1,:)-XY(2,:)); % length of line
    dx=XY(2,1)-XY(1,1);
    dy=XY(2,2)-XY(1,2);
    alpha=atan2(dy,dx); % angle of rotation
    cosa=cos(alpha);
    sina=sin(alpha);

    points=linspace(pi/4,3*pi/4);
    a=(0.5-cos(points)/2^0.5)*l;
    b=((sin(points)-2^0.5/2)/(1-2^0.5/2))*Radius;

    x=a*cosa-b*sina+XY(1,1);
    y=a*sina+b*cosa+XY(1,2);

    % Arrow end
    xai=x([end end-20]);
    yai=y([end end-20]);
    xa1=LengthArrow*[0 0.1 0.08 0.1 0]';
    ya1=LengthArrow*[0 0.03 0 -0.03 0]';
    dx=diff(xai);
    dy=diff(yai);
    alpha=atan2(dy,dx); % angle of rotation
    cosa=cos(alpha);
    sina=sin(alpha);
    xa=xa1*cosa-ya1*sina+xai(1);
    ya=xa1*sina+ya1*cosa+yai(1);

    plot(x,y,'-','color',Color,'linewidth',LineWidth); hold on
    fill(xa,ya,repmat(Color,1,1))
    plot(xa,ya,'-','color',Color)

end




%% Update
function Update(handles)
    tic
    global green_image step_execution gray_image waiting_image
    % Get data
    Processing('msj','Getting Data',handles)
    Experiment=Get_Experiment(handles);
    if isfield(Experiment,'Data')
        if isfield(Experiment.Data,'Data')
            set(handles.PlotDataPushbutton,'CData', gray_image);
            set(handles.PlotPeaksPushbutton,'CData', gray_image);
            set(handles.DimRedPushbutton,'CData', gray_image);
            set(handles.ClusteringPushbutton,'CData', gray_image);
            set(handles.PlotNetworkPushbutton,'CData', gray_image);
            step_execution = [0 0 0 0 0];
            
            set(handles.PlotDataPushbutton,'CData', waiting_image);
            
            Processing('msj','Plotting Simple Experiment',handles)
            Experiment=PlotSimpleExperiment(Experiment,handles);
            Processing('end','Plotting Simple Experiment',handles)
            
            
            Processing('msj','Plotting Activity',handles)
            Experiment=PlotActivity(Experiment,handles);
            Processing('end','Plotting Activity',handles)
            
            set(handles.PlotDataPushbutton,'CData', green_image);
            step_execution = [1 0 0 0 0];
            
            set(handles.PlotPeaksPushbutton,'CData', waiting_image);
            
            Processing('msj','Plotting Peaks',handles)
            Experiment=PlotPeaks(Experiment,handles);
            Processing('end','Plotting Peaks',handles)
            
            set(handles.PlotPeaksPushbutton,'CData', green_image);
            step_execution = [1 1 0 0 0];
            
            set(handles.DimRedPushbutton,'CData', waiting_image);
            
            Processing('msj','Plotting Dimensional Reduction',handles)
            Experiment=PlotDimRed(Experiment,handles);
            Processing('end','Plotting Dimensional Reduction',handles)
            
            set(handles.DimRedPushbutton,'CData', green_image);
            step_execution = [1 1 1 0 0];
            
            set(handles.ClusteringPushbutton,'CData', waiting_image);
            
            Processing('msj','Plotting Clustering',handles)
            Experiment=PlotClustering(Experiment,handles);
            Processing('end','Plotting Clustering',handles)
            Processing('msj','Plotting LLE',handles)
            Experiment=PlotHtoLLE(Experiment,handles);
            
            set(handles.ClusteringPushbutton,'CData', green_image);
            step_execution = [1 1 1 1 0];
            
            set(handles.PlotNetworkPushbutton,'CData', waiting_image);
            
            Processing('msj','Plotting Group Structure',handles)
            Experiment=PlotGroupStruc(Experiment,handles);
            Processing('end','Plotting Group Structure',handles)
            
            Processing('msj','Plotting Adjacency Matrix',handles)
            Experiment=PlotAdjacent(Experiment,handles);
            Processing('end','Plotting Adjacency Matrix',handles)
            
            Processing('msj','Plotting Network',handles)
            Experiment=PlotNetwork(Experiment,handles);
            Processing('end','Plotting Adjacency Matrix',handles)
            
             set(handles.PlotNetworkPushbutton,'CData', green_image);
             step_execution = [1 1 1 1 1];
            
            Processing('end','Analysis Completed',handles)
            assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
        else
            Processing('error','Getting Structs',handles)
        end
    else
        Processing('error','Getting Data',handles)
    end
    toc
end
%% ...::: Plot functions :::...
%% Plot Simple
function Experiment=PlotSimpleExperiment(Experiment, handles)
    
    % Inputs
    Data=Experiment.Data.Data;
    DataName=Experiment.Data.DataName;
    XYColors=Experiment.Colors.XY;
    F=Experiment.Data.Frames;
    C=Experiment.Data.Neurons;
   
    % Plot experiment data as a raster, plot activity and coactivity.
    
    % Figure colors 
    BGColor=[1 1 1];
    ColorLight=[0.7 0.7 0.7];
    ColorDark=[0 0 0];
    ColorTextLines=[0 0 0];

    % Plot data
    obj=findobj('name',['Experiment (' DataName ')']);
    if isempty(obj)
        figure('name',['Experiment (' DataName ')'],'numbertitle','off','position',[0 0 1220 460]);
    else
        figure(obj);
    end
    clf
    axes('outerposition',[0.125 0.34 0.75 0.66])
    set(gca,'Tag',['RasterAxes' DataName])
    axis([1 F 1 C]); box; hold on
    for i=1:C
        idx=find(Data(:,i));
        plot(idx,Data(idx,i)*i,'.','color',XYColors(i,:))
    end
    hold off
    title(DataName)
    %set(gca,'color',BGColor,'XTicklabel','','XTick',0,'YTick',0,...
    %set(gca,'color',BGColor,'XTicklabel','','XTick',[1 C],'YTick',[1 C],...
    %    'xcolor',ColorDark,'ycolor',ColorDark)
    set(gca,'color',BGColor,'XTicklabel','','XTick',[1 C],...
        'xcolor',ColorDark,'ycolor',ColorDark)
    
    % Plot Coactive neurons
    axes('outerposition',[0.125 0 0.75 0.34]);
    Co=sum(Data,2);
    fpm=str2num(get(handles.FPSEdit,'string'))*60;
    %{
    [S,Fq,T,P] = spectrogram(Co,50,49,512,max(Co)*2); % max(Co)*2 for fitting plot 
    surf(T*max(Co)*2/fpm,Fq,10*log10(P),'edgecolor','none'); axis tight; hold on;
    view(0,90);
    %}
    
    plot((1:F)/fpm,Co,'k')
    title('Coactive neurons')
    set(gca,'Tag',['CoactiveAxes' DataName])
    set(gca,'color',BGColor,'box','off','xcolor',ColorDark,'ycolor',ColorDark)
    set(gca,'xtick',0:(F/fpm))
    xlabel('Time (min)'); xlim([0 F/fpm])
    ylabel('# active neurons'); ylim([0 max(Co)])
    
    % Plot Accumulated activity
    axes('outerposition',[0.875 0.34 0.125 0.66]);
    Ac=sum(Data,1)*100/F;
    plot(Ac,'-k')
    hold on
   
    % Settings plot
    title('Activity')
    set(gca,'XTicklabel','','XTick',0,'xlim',[0 C],'view',[90 -90])
    set(gca,'Tag',['ActivityAxes' DataName])
    set(gca,'color',BGColor,'box','off','XTick',0,...
    'xcolor',ColorDark,'ycolor',ColorDark)
    ylabel('% Active frames'); ylim([0 max(Ac)])
    
    % Output
    Experiment.Raster.Activity=Ac;
    Experiment.Raster.Coactivity=Co; 
end

%% Plot Activity 
function Experiment=PlotActivity(Experiment,handles)

    Activity=Experiment.Raster.Activity;
    DataName=Experiment.Data.DataName;
    
    obj=findobj('name',['Activity (' DataName ')']);
    if isempty(obj)
        figure('name',['Activity (' DataName ')'],'numbertitle','off','position',[0 0 900 300]);
    else
        figure(obj);
    end
    clf
    
    % Figure colors 
    BGColor=[1 1 1];
    ColorLight=[0.9 0 0.7];
    ColorDark=[0.5 0 0.3];
    ColorTextLines=[0 0 0];
    
    % Plot activity histogram
    axes('outerposition',[0 0 0.5 1])
    set(gca,'Tag',['HistAxes' DataName])
    N=length(Activity);
    nbins=round(log2(N)+1);
    [y x]=hist(Activity,nbins);
    idx=find(y);
    y=y(idx);
    x=x(idx);
    %bar(x,y,'facecolor',ColorLight,'edgecolor',ColorDark)
    plot(x,y,'.','color',ColorLight,'MarkerSize',30)
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    hold on
    
    % Plot power law fit
    [slope, intercept, R2] = logfit(x,y,'loglog');
    xfit=min(x):(max(x)-min(x))/100:max(x);
    yfit=(10^intercept)*xfit.^(slope);
    plot(xfit,yfit,'--','color',ColorDark,'linewidth',2.5)
    
    % Plot labels and colors
    title('Activity histogram','color',ColorTextLines)
    xlabel('% Active frames','color',ColorTextLines)
    ylabel('Count of neurons','color',ColorTextLines)
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    
    % Write the slope of power law fitted (Power law:
    % y=10^intercept*x^slope)
    text(max(x)/2,max(y)*0.8,['\alpha=' num2str(-slope,'%0.1f')],'color',ColorDark,'FontSize',14,'FontWeight','bold')
    text(max(x)/2,max(y)*0.6,['R^2=' num2str(R2,'%0.3f')],'color',ColorDark,'FontSize',14,'FontWeight','bold')
    
    % Loglog activity histogram
    axes('outerposition',[0.5 0 0.5 1])
    loglog(x,y,'.','color',ColorLight,'MarkerSize',30)
    hold on
    loglog(xfit,yfit,'--','color',ColorDark,'linewidth',2.5)
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    
    % Plot labels and colors
    title('Log-Log Activity histogram','color',ColorTextLines)
    xlabel('% Active frames','color',ColorTextLines)
    ylabel('','color',ColorTextLines)
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    set(gca,'Tag',['LogAxes' DataName])
    
    % Outputs
    Experiment.Network.ActivitySlope=-slope;
    
end

%% Plot Peaks
function Experiment=PlotPeaks(Experiment,handles)
    
    % Inputs
    Data=Experiment.Data.Data;
    DataName=Experiment.Data.DataName;
    Co=Experiment.Raster.Coactivity;
    Ac=Experiment.Raster.Activity;
    
    % Colors
    BGColor=[1 1 1];
    ColorLight=[0.5 0 0];
    ColorDark=[0.3 0 0];
    ColorTextLines=[0 0 0];
    
    % Get Peaks
    ThSync=str2num(get(handles.PeaksThEdit,'string'));
    
    % Find peaks
    Join=get(handles.JoinPeaksCheckbox,'value');
    PeaksIdx=FindPeaks_JP(Co,ThSync,Join);
    
    PeaksIdx_Join=FindPeaks_JP(Co,ThSync,true);
    idx=find(PeaksIdx_Join);
    PeaksIdx_Join=PeaksIdx_Join(idx);
    
    Binary=get(handles.BinaryPeaksCheckbox,'value');
    DataPeaks=PeaksVectors_JP(Data,PeaksIdx,~Binary);
    [ActiveDataPeaks NoActive Active]=RActive_JP(DataPeaks,0);
    
    % Similarity
    Num=get(handles.SimMethodPopupmenu,'Value');
    Strings=get(handles.SimMethodPopupmenu,'String');
    SimMethod=Strings{Num};
    switch SimMethod
        case 'LLEdistance'
            X=ActiveDataPeaks';
            N=size(X,2);
            X2=sum(X.^2,1);
            Distance=repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
        otherwise
            Distance=squareform(pdist(ActiveDataPeaks,SimMethod));
    end
    Sim=1-Distance/max(max(Distance));
    
    % Figure Peaks
    obj=findobj('name',['Peaks (' DataName ')']);
    if isempty(obj)
        figure('name',['Peaks (' DataName ')'],'numbertitle','off','position',[0 0 300 300]);
    else
        figure(obj);
    end
    clf

    % Plot similarity
    imagesc(Sim)
    set(gca,'YDir','normal','xcolor',ColorDark,'ycolor',ColorDark)
    title([SimMethod ' similarity'],'color',ColorDark)
    xlabel('# peak (t)')
    set(gca,'Tag',['SimilarityAxes' DataName])
    
    % Indexes of sort Activity
    [value ActivityRank]=sort(Ac(Active),'descend');
    
    % Outputs
    Experiment.Peaks.Threshold=ThSync;
    Experiment.Peaks.Join=Join;
    Experiment.Peaks.Binary=Binary;
    Experiment.Peaks.Method=SimMethod;
    Experiment.Peaks.Similarity=Sim;
    Experiment.Peaks.ActiveNeurons=Active;
    Experiment.Peaks.ActiveDataPeaks=ActiveDataPeaks;
    Experiment.Peaks.Index=PeaksIdx;
    Experiment.Ranks.Activity=ActivityRank;
    Experiment.Peaks.Join=PeaksIdx_Join;
end

%% Plot Dimensional Reduction
function Experiment=PlotDimRed(Experiment,handles)
    
    % Inputs
    DataName=Experiment.Data.DataName;
    ActiveDataPeaks=Experiment.Peaks.ActiveDataPeaks;
    Distance=1-Experiment.Peaks.Similarity;
    SimMethod=Experiment.Peaks.Method;
                
    RED=get(handles.RedRadiobutton,'value');
    if RED
        % Figure Dimensional Reduction
        obj=findobj('name',['Dimensional Reduction (' DataName ')']);
        if isempty(obj)
            figure('name',['Dimensional Reduction (' DataName ')'],'numbertitle','off','position',[0 0 600 300]);
        else
            figure(obj);
        end
        clf
        
        % Axes positions
        Pos1=[0 0 0.5 1];
        Pos2=[0.5 0 0.5 1];
        
        % Colors
        BGColor=[1 1 1];
        ColorLight=[0.5 0 0];
        ColorDark=[0.3 0 0];
        ColorTextLines=[0 0 0];
        
        % Reduction
        Num=get(handles.DimRedPopupmenu,'Value');
        Strings=get(handles.DimRedPopupmenu,'String');
        Reduction=Strings{Num};
        switch Reduction
            case 'LLE'
                Neighbors=get(handles.NeighborsEdit,'String');
                Neighbors=str2num(Neighbors);
                if Neighbors<=size(ActiveDataPeaks,1)
                    DataRed=LLE_JP(ActiveDataPeaks',Neighbors,3,Distance);
                else
                    Neighbors=1;
                    set(handles.NeighborsEdit,'String','1');
                    DataRed=LLE_JP(ActiveDataPeaks',Neighbors,3,Distance);
                end
            case 'PCA'
                DataRed=princomp(ActiveDataPeaks');
                DataRed=DataRed(:,1:3)';
                Neighbors=[];
        end
        
        % Plot Data Reduced
        axes('outerposition',Pos1)
        plot3(DataRed(1,:),DataRed(2,:),DataRed(3,:),'.','color',ColorDark,'MarkerSize',20)
        xl=get(gca,'xlim'); yl=get(gca,'ylim'); zl=get(gca,'zlim');
        delta=max([abs(xl(1)-xl(2)) abs(yl(1)-yl(2)) abs(zl(1)-zl(2))]);
        xl(2)=xl(1)+delta; yl(2)=yl(1)+delta; zl(2)=zl(1)+delta;
        xlabel([Reduction ' 1']);ylabel([Reduction ' 2']);zlabel([Reduction ' 3']);
        set(gca,'xlim',xl,'ylim',yl,'zlim',zl)
        set(gca,'XTick',[],'YTick',[],'ZTick',[])
        set(gca,'box','on','view',[0 90])
        title([Reduction ' reduction'])
        set(gca,'Tag',['DimRedAxes1' DataName])
        
        % Similarity
        switch SimMethod
        case 'LLEdistance'
            X=DataRed;
            N=size(X,2);
            X2=sum(X.^2,1);
            DistanceR=repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;
        otherwise
            DistanceR=squareform(pdist(DataRed',SimMethod));
        end
        SimR=1-DistanceR/max(max(DistanceR));
        
        % Plot similarity
        axes('outerposition',Pos2)
        imagesc(SimR)
        set(gca,'YDir','normal','xcolor',ColorDark,'ycolor',ColorDark)
        title([SimMethod ' similarity (data reduced)'],'color',ColorDark)
        xlabel('# peak (t)')
        set(gca,'Tag',['SimilarityRedAxes' DataName])
        
    else
        Reduction='No reduction';
        Neighbors=[];
        DataRed=[];
        SimR=[];
    end
    
    % Outputs
    Experiment.DimensionalReduction.Method=Reduction;
    Experiment.DimensionalReduction.Neighbors=Neighbors;
    Experiment.DimensionalReduction.ActiveDataPeaks=DataRed;
    Experiment.DimensionalReduction.Similarity=SimR;
end

%% Plot Clustering
function Experiment=PlotClustering(Experiment,handles)
    
    % Inputs
    DataName=Experiment.Data.DataName;
    RED=get(handles.RedRadiobutton,'value');
    HCA=get(handles.HierarchicalRadiobutton,'value');
    Groups=str2num(get(handles.GroupsThEdit,'string'));
    CleanTh=str2num(get(handles.CleanThEdit,'string'));
    Sim=Experiment.Peaks.Similarity;
    SimRed=Experiment.DimensionalReduction.Similarity;
    SimMethod=Experiment.Peaks.Method;
    DataRed=Experiment.DimensionalReduction.ActiveDataPeaks;
    GroupColors=Experiment.Colors.States;
    Reduction=Experiment.DimensionalReduction.Method;
    
    
    % Dimensional reduction
    if RED
        Sim2Clust=SimRed;
    else
        Sim2Clust=Sim;
    end
    
    % Clustering method
    warning off
    if HCA
        % Clustering linkage method
        ClustMethod='Hierarchical';
        Num=get(handles.ClustPopupmenu,'Value');
        Strings=get(handles.ClustPopupmenu,'String');
        LinkageMethod=Strings{Num};
        Tree=linkage(squareform(1-Sim2Clust,'tovector'),LinkageMethod);
        idx_vector_state=cluster(Tree,'maxclust',Groups);
        Experiment.Clustering.Tree=Tree;
    else
        ClustMethod='K-means';
        LinkageMethod='None';
        ActiveDataPeaks=Experiment.Peaks.ActiveDataPeaks;
        idx_vector_state=kmeans(ActiveDataPeaks,Groups,'emptyaction','singleton');
    end
    warning on

    % Figure Clustering
    % Axes positions
    if HCA
        if RED
            size_figure=[0 0 1200 600];
            PosSim=[0 0.5 0.25 0.5]; PosSimRed=[0.25 0.5 0.25 0.5];
            PosSort=[0.5 0.5 0.25 0.5]; PosTree=[0.75 0.5 .25 0.5];
            PosData=[0.375 0 0.25 0.5];
        else
            size_figure=[0 0 900 600];
            PosSim=[0 0.5 0.33 0.5]; PosSort=[0.33 0.5 0.33 0.5];
            PosTree=[0.66 0.5 0.33 0.5]; PosData=[0 0 .33 0.5];
        end
    else
        if RED
            size_figure=[0 0 900 600];
            PosSim=[0 0.5 0.33 0.5]; PosSimRed=[0.33 0.5 0.33 0.5];
            PosSort=[0.66 0.5 0.34 0.5]; PosData=[0.33 0 .34 0.5];
        else
            size_figure=[0 0 900 300];
            PosSim=[0 0 0.33 1]; PosSort=[0.33 0 0.33 1];
            PosData=[0.66 0 0.34 1];
        end
    end
    
    obj=findobj('name',['Clustering (' DataName ')']);
    if isempty(obj)
        figure('name',['Clustering (' DataName ')'],'numbertitle','off','position',size_figure);
    else
        figure(obj);
    end
    clf
        
    % Colors
    BGColor=[1 1 1];
    ColorLight=[0.5 0 0];
    ColorDark=[0.3 0 0];
    ColorTextLines=[0 0 0];
    
    % Plot Similarity from original data
    axes('outerposition',PosSim)
    imagesc(Sim)
    set(gca,'YDir','normal','xcolor',ColorDark,'ycolor',ColorDark)
    title([SimMethod ' similarity'],'color',ColorDark)
    xlabel('# peak (t)')
    ratio=[];
    
    T=[];
    PeaksIdx=Experiment.Peaks.Index;
    Peaks=Join_Peaks(idx_vector_state,PeaksIdx);
    
    if RED
        % Plot Similarity from original data
        axes('outerposition',PosSimRed)
        imagesc(SimRed)
        set(gca,'YDir','normal','xcolor',ColorDark,'ycolor',ColorDark)
        title([SimMethod ' similarity  (data reduced)'],'color',ColorDark)
        xlabel('# peak (t)')
        
        % Plot Data
        axes('outerposition',PosData)
        for i=1:Groups
            idx=find(idx_vector_state==i);
            plot3(DataRed(1,idx),DataRed(2,idx),DataRed(3,idx),'.','color',GroupColors(i,:),'MarkerSize',20)
            hold on
        end
        xl=get(gca,'xlim'); yl=get(gca,'ylim'); zl=get(gca,'zlim');
        delta=max([abs(xl(1)-xl(2)) abs(yl(1)-yl(2)) abs(zl(1)-zl(2))]);
        xl(2)=xl(1)+delta; yl(2)=yl(1)+delta; zl(2)=zl(1)+delta;
        xlabel([Reduction ' 1']);ylabel([Reduction ' 2']);zlabel([Reduction ' 3']);
        set(gca,'xlim',xl,'ylim',yl,'zlim',zl)
        set(gca,'XTick',[],'YTick',[],'ZTick',[])
        set(gca,'box','on','view',[0 90])
        title([Reduction ' reduction'])
        set(gca,'Tag',['DimRedAxes1' DataName])
    else
        % Plot Data
        axes('outerposition',PosData)
        PeaksIdx=Experiment.Peaks.Index;
        Peaks=Join_Peaks(idx_vector_state,PeaksIdx);
        XY=[cos(2*pi*[1:Groups]'/Groups) sin(2*pi*[1:Groups]'/Groups)];
        for i=1:Groups
            idx=find(idx_vector_state==i);
            idx2=find(Peaks==i);
            ratio(i)=length(idx2)/length(Peaks)*100;
            sd_g=std(1-squareform(1-Sim(idx,idx),'tovector'));
            mean_g=mean(1-squareform(1-Sim(idx,idx),'tovector'));
            cv_g=sd_g/mean_g;
            if isnan(cv_g)
                cv_g=0;
            end
            %plot(XY(i,1),XY(i,2),'.','color',[0.3 0.3 0.3],'MarkerSize',100*cv_g+30);
            plot(XY(i,1),XY(i,2),'.','color',GroupColors(i,:),'MarkerSize',100*cv_g+80)
            hold on
            text(XY(i,1),XY(i,2),[num2str(ratio(i),'%2.0f') '%'],...
                'HorizontalAlignment','Center','FontWeight','bold','FontSize',12)
        end
        set(gca,'box','off','xtick',[],'ytick',[],'xcolor',BGColor,...
            'ycolor',BGColor,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5])
        title('Ratio of states')
        
        % Plot transitions
        T=StatesTrans(Peaks);
        if Groups==3
            if HCA
                Radius=-0.15;
                LengthArrow=3.5;
                LineWidth=2;
                ColorArrow=[0.5 0.5 0.5];
                TextColor=[0 0 0];
                
                % Ratio of all transitions
                axes('outerposition',[0.33 0 0.33 .5])
                % Loops
                for i=1:3
                    plot(XY(i,1),XY(i,2),'.','color',GroupColors(i,:),'MarkerSize',100)
                    hold on
                    if T(i,i)
                        text(XY(i,1)*1.4,XY(i,2)*1.4,[num2str(floor(T(i,i)*100)) '%'],'HorizontalAlignment','center',...
                        'FontSize',14,'FontWeight','bold','color',TextColor)
                    end
                end
                
                % Arrows
                Trans=[1 2;2 3;3 1;2 1;3 2;1 3];
                Alignment={'right','left','left','center','center','center'};
                for i=1:6
                    from=Trans(i,1);
                    to=Trans(i,2);
                    if T(from,to)
                        [x y]=PlotArrow(XY([from to],:),Radius,LengthArrow,LineWidth,ColorArrow);
                        text(x(50),y(50),[num2str(floor(T(from,to)*100)) '%'],...
                            'HorizontalAlignment',Alignment{i},...
                            'FontSize',14,'FontWeight','bold','color',TextColor)
                    end
                end
                set(gca,'box','off','xtick',[],'ytick',[],'xcolor',BGColor,...
                    'ycolor',BGColor,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5])
                %title('Ratio of all transitions')
                Entropy=0;
                for i=1:9
                    if T(i)
                        Entropy=Entropy+T(i)*log2(T(i));
                    end
                end
                title(['Ratio of all transitions; Entropy=' num2str(-Entropy) ' bits'])
                
                % Ratio of state transitions
                for i=1:3
                    if sum(T(i,:))
                        T(i,:)=T(i,:)./sum(T(i,:));
                    end
                end
                axes('outerposition',[0.66 0 0.33 0.5])
                % Loops
                for i=1:3
                    plot(XY(i,1),XY(i,2),'.','color',GroupColors(i,:),'MarkerSize',100)
                    hold on
                    if T(i,i)
                        text(XY(i,1)*1.4,XY(i,2)*1.4,[num2str(floor(T(i,i)*100)) '%'],'HorizontalAlignment','center',...
                        'FontSize',14,'FontWeight','bold','color',TextColor)
                    end
                end
                
                % Arrows
                Trans=[1 2;2 3;3 1;2 1;3 2;1 3];
                Alignment={'right','left','left','center','center','center'};
                for i=1:6
                    from=Trans(i,1);
                    to=Trans(i,2);
                    if T(from,to)
                        [x y]=PlotArrow(XY([from to],:),Radius,LengthArrow,LineWidth,ColorArrow);
                        text(x(50),y(50),[num2str(floor(T(from,to)*100)) '%'],...
                            'HorizontalAlignment',Alignment{i},...
                            'FontSize',14,'FontWeight','bold','color',TextColor)
                    end
                end
                set(gca,'box','off','xtick',[],'ytick',[],'xcolor',BGColor,...
                    'ycolor',BGColor,'xlim',[-1.5 1.5],'ylim',[-1.5 1.5])
                title('Ratio of state transitions')
                
            end
        end
    end
    
    if HCA
        % Plot dendrogram
        axes('outerposition',PosTree)
        dendrogram(Tree,0,'orientation','rigth');%,'colorthreshold',2); %2 is for figure
        Tid=str2num(get(gca,'YTicklabel'));
        set(gca,'xtick',[],'ytick',[],'xcolor',BGColor,'ycolor',BGColor)
        title([LinkageMethod ' linkage'])
        set(gca,'Tag',['DendroAxes' DataName])

        % Plot sort by similarity
        axes('outerposition',PosSort)
        imagesc(Sim2Clust(Tid,Tid))
        set(gca,'YDir','normal','xcolor',ColorDark,'ycolor',ColorDark)
        %set(gca,'YTicklabel',Tid,'XTicklabel',Tid,...
        %    'YTick',[1:length(Sim2Clust)],'XTick',[1:length(Sim2Clust)])
        set(gca,'xtick',[],'ytick',[])
        title('Sort by similarity')
        set(gca,'Tag',['SimSortAxes' DataName])
    else
        % Sort vector index by states
        Tid=[];
        for i=1:Groups
            idx=find(idx_vector_state==i);
            Tid=[Tid; idx];
        end
        
        % Plot sort by similarity
        axes('outerposition',PosSort)
        imagesc(Sim2Clust(Tid,Tid))
        set(gca,'YDir','normal','xcolor',ColorDark,'ycolor',ColorDark)
        %set(gca,'YTicklabel',Tid,'XTicklabel',Tid,...
        %    'YTick',[1:length(Sim2Clust)],'XTick',[1:length(Sim2Clust)])
        set(gca,'xtick',[],'ytick',[])
        title('Sort by similarity')
        set(gca,'Tag',['SimSortAxes' DataName])
    end
    
    % Plot LLE
    %Join=NMDA_Analysis.Peaks.Join;
    
    % Outputs
    Experiment.Clustering.Method=ClustMethod;
    Experiment.Clustering.LinkageMethod=LinkageMethod;
    Experiment.Clustering.TotalStates=Groups;
    Experiment.Clustering.VectorStateIndex=idx_vector_state;
    Experiment.Clustering.StatesRatio=ratio;
    Experiment.Clustering.StatesSequence=Peaks;
    Experiment.Clustering.TotalPeaks=length(Peaks);
    Experiment.Clustering.TransitionsRatio=T;
    
end

%% Plot Hierarchical to LLE
function Experiment=PlotHtoLLE(Experiment,handles)
    
    % Inputs
    DataName=Experiment.Data.DataName;
    States=Experiment.Clustering.VectorStateIndex;
    Peaks=Experiment.Peaks.ActiveDataPeaks;
    Join=Experiment.Peaks.Join;
    GroupColors=Experiment.Colors.States;
    
    % Join Peaks
    PeaksJoin=[];
    if (size(Peaks,1)==size(Join,1))
        for i=1:max(Join)
            idx=find(Join==i);
            PeaksJoin(i,:)=sum(Peaks(idx,:),1);
            StatesJoin(i)=States(idx(1));
        end
    else
        PeaksJoin=Peaks;
        StatesJoin=States;
    end
    Groups=max(StatesJoin);

    % LLE distance
    X=PeaksJoin';
    N=size(X,2);
    X2=sum(X.^2,1);
    Distance=repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

    if N>50
        Max_N=50;
    else
        Max_N=N-1;
    end

    % Figure Dimensional Reduction
    obj=findobj('name',['Optimal Neighbors (' DataName ')']);
    if isempty(obj)
        figure('name',['Optimal Neighbors (' DataName ')'],'numbertitle','off','position',[0 0 600 300]);
    else
        figure(obj);
    end
    clf

    warning off
    DIs=[];
    for Neighbors=2:Max_N
        try
            DataRed=LLE_JP(PeaksJoin',Neighbors,3,Distance);
            DistanceRed=squareform(pdist(DataRed','Euclidean'));
            DI = DunnIdx_JP(Groups,DistanceRed,StatesJoin);
            DIs(Neighbors)=DI;
        catch ME
            disp(['No LLE for ' num2str(Neighbors) ' neighbors.'])
        end
    end

    plot(DIs); hold on;
    [val NeighborsIDEAL]=max(DIs);
    plot(NeighborsIDEAL,val,'*r','MarkerSize',20)
    try
        DataRed=LLE_JP(PeaksJoin',NeighborsIDEAL,3,Distance);
        DistanceRed=squareform(pdist(DataRed','Euclidean'));
        title(['Optimal Neighbors: ' num2str(NeighborsIDEAL) ])
        warning on

        % Figure Dimensional Reduction
        obj=findobj('name',['Dimensional Reduction (' DataName ')']);
        if isempty(obj)
            figure('name',['Dimensional Reduction (' DataName ')'],'numbertitle','off','position',[0 0 600 300]);
        else
            figure(obj);
        end
        clf

        % Axes positions
        Pos1=[0 0 0.5 1];
        Pos2=[0.5 0 0.5 1];

        % Colors
        BGColor=[1 1 1];
        ColorLight=[0.5 0 0];
        ColorDark=[0.3 0 0];
        ColorTextLines=[0 0 0];

        % Plot Data Reduced
        axes('outerposition',Pos1)
        plot3(DataRed(1,:),DataRed(2,:),DataRed(3,:),'.','color',ColorDark,'MarkerSize',20)
        xl=get(gca,'xlim'); yl=get(gca,'ylim'); zl=get(gca,'zlim');
        delta=max([abs(xl(1)-xl(2)) abs(yl(1)-yl(2)) abs(zl(1)-zl(2))]);
        xl(2)=xl(1)+delta; yl(2)=yl(1)+delta; zl(2)=zl(1)+delta;
        xlabel('LLE 1');ylabel('LLE 2');zlabel('LLE 3');
        set(gca,'xlim',xl,'ylim',yl,'zlim',zl)
        set(gca,'XTick',[],'YTick',[],'ZTick',[])
        set(gca,'box','on','view',[0 90])
        title('LLE reduction')
        set(gca,'Tag',['DimRedAxes1' DataName])

        % Plot Data
        axes('outerposition',Pos2)
        for i=1:Groups
            idx=find(StatesJoin==i);
            plot3(DataRed(1,idx),DataRed(2,idx),DataRed(3,idx),'.','color',GroupColors(i,:),'MarkerSize',40)
            hold on
        end
        xlabel('LLE 1');ylabel('LLE 2');zlabel('LLE 3');
        
        % Outputs
        Experiment.Clustering.HtoLLE=DataRed;
    end
end

%% Plot Groups Structure
function Experiment=PlotGroupStruc(Experiment,handles)
    
    % Inputs
    DataName=Experiment.Data.DataName;
    CleanTh=str2num(get(handles.CleanThEdit,'string'));
    Groups=Experiment.Clustering.TotalStates;
    GroupColors=Experiment.Colors.States;
    PeaksIdx=Experiment.Peaks.Index;
    Active=Experiment.Peaks.ActiveNeurons;

    % Get states of each neuron
    Experiment=GetStatesNeuron(Experiment);
    CellsGroups=Experiment.Structure.NeuronsStates;
    CellsActGroups=Experiment.Structure.RateNeuronsStates;
    
    % Clean outliers
    C=length(Active);
    tot=sum(CellsActGroups,1);
    CellsActGroupsNormClean=CellsActGroups;
    for i=1:Groups
        Norms=find(CellsActGroups(i,:)./tot*100>=CleanTh);
        for j=1:length(Norms)
            CellsActGroupsNormClean(:,Norms(j))=0;
            CellsActGroupsNormClean(i,Norms(j))=1;
        end
    end
    CellsGroupsClean=CellsActGroupsNormClean;
    for i=1:Groups
        idx=find(CellsActGroupsNormClean(i,:));
        CellsGroupsClean(i,idx)=i+1;
    end
    Experiment.Structure.NeuronsStates=CellsGroupsClean;
    
    % Get states combinations of each neuron
    Experiment=GetStatesComb(Experiment);
    StatesClean=Experiment.Structure.StatesRatio;
    LabelsClean=Experiment.Structure.LabelsStates;
    IndexesClean=Experiment.Structure.Indexes;
    IndexesCombClean=Experiment.Structure.StatesCombinations;
    idxSortClean=Experiment.Structure.IndexesLabelsSort;
    
    % Plot Groups structure
    obj=findobj('name',['Groups structure (' DataName ')']);
    if isempty(obj)
        figure('name',['Groups structure (' DataName ')'],'numbertitle','off','position',[0 0 1200 600]);
    else
        figure(obj);
    end
    clf
    
    % Colors
    BGColor=[1 1 1];
    ColorLight=[0.9 0 0.7];
    ColorDark=[0.5 0 0.3];
    ColorTextLines=[0 0 0];
    
    % Axes Posistion
    aGroupsClean=[0 0.8 0.75 0.2];
    aGroupsCleanSort=[0 0.6 0.75 0.2];
    aPeaksClean=[0 0.4 0.75 0.2];
    aRatioStatesClean=[0.75 0.45 0.25 0.5];
    
    aHistPeaks=[0.125 0 0.375 0.4];
    aHistPeaksLog=[0.5 0 0.375 0.4];
    
    % 1. Plot # neuron VS States
    axes('outerposition',aGroupsClean)
    image(CellsGroupsClean(:,IndexesClean))
    colormap([1 1 1; GroupColors])
    set(gca,'Tag',['GroupAxes' DataName])
    set(gca,'YDir','normal','xcolor',ColorTextLines,'ycolor',ColorTextLines,'YTick',1:Groups)
    set(gca,'Xtick',[])
    ylabel('State')
    title('Groups structure')
    hold on
    for i=1:Groups-1
        plot([0 C+2],[i+.5 i+.5],'-k')
    end
    for i=1:C-1
        plot([i+0.5 i+0.5],[0 Groups+0.5],'-k')
    end
    
    % 2. Plot # neuron VS % of each state
    axes('outerposition',aGroupsCleanSort);
    tot=sum(CellsActGroupsNormClean,1);
    for j=Groups:-1:1
        bar(1:C,sum(CellsActGroupsNormClean(1:j,IndexesClean),1)./tot(IndexesClean)*100,...
            'facecolor',GroupColors(j,:),'edgecolor',[0 0 0])
        hold on
    end
    set(gca,'Xlim',[.5 C+.5],'XTick',[])
    ylabel('% in each state')
    
    % 3. Plot # neuron VS % of each state
    axes('outerposition',aPeaksClean);
    hold on
    Comb=2^Groups-1;
    for i=1:Comb
        bin=dec2bin(i,Groups);
        color=[];
        for j=1:Groups
            if str2num(bin(j))
                color=[color; GroupColors(j,:)];
            end
        end
        ColorComb(i,:)=mean(color,1);
    end
    Peaks=max(PeaksIdx);
    maxCellPeaks=max(sum(CellsActGroups,1))/Peaks*100;
    plot([0 C],[maxCellPeaks*.3 maxCellPeaks*.3],'--k')
    for i=1:C
        bar(i,sum(CellsActGroups(:,IndexesClean(i)),1)/Peaks*100,...
            'facecolor',ColorComb(IndexesCombClean(IndexesClean(i)),:),'edgecolor',[0 0 0])
    end
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    set(gca,'Xlim',[.5 C+.5],'XTick',1:C,'XTickLabel',Active(IndexesClean))
    set(gca,'Ylim',[0 maxCellPeaks])
    xlabel('# neuron');ylabel('% peaks')
    [value Experiment.Ranks.Peaks]=sort(sum(CellsActGroups,1),'descend');
    
    % 4. Plot % of states combinations (after clean)
    axes('outerposition',aRatioStatesClean); hold on
    for i=1:Comb
        bar(i,StatesClean(i)*100,'facecolor',ColorComb(idxSortClean(i),:),'edgecolor',ColorTextLines)
        hold on
        text(i,StatesClean(i)*100+1,LabelsClean{i},'HorizontalAlignment','Center')
    end
%     bar(StatesClean1*100,'facecolor',ColorLight,'edgecolor',ColorDark)
%         text([1:length(StatesClean1)],StatesClean1*100+1,LabelsClean1,'HorizontalAlignment','Center')
    xlabel('States','color',ColorTextLines)
    ylabel('% of neurons','color',ColorTextLines)
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    set(gca,'XTick',[])
    
    % 5. Plot peaks histogram
    axes('outerposition',aHistPeaks);
    PeaksVSNeuron=sum(CellsActGroups(:,IndexesClean),1);
    N=length(PeaksVSNeuron);
    nbins=round(log2(N)+1);
    [y x]=hist(PeaksVSNeuron,nbins);
    idx=find(y);
    y=y(idx);
    x=x(idx);
    plot(x,y,'.','color',ColorLight,'MarkerSize',30); hold on
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    xlabel('# peaks')
    ylabel('Count of neurons')
    %title('Peaks histogram','color',ColorTextLines)
    % Plot power law fit
    [slope, intercept, R2] = logfit(x,y,'loglog');
    xfit=min(x):(max(x)-min(x))/100:max(x);
    yfit=(10^intercept)*xfit.^(slope);
    plot(xfit,yfit,'--','color',ColorDark,'linewidth',2.5)
    % Write the function fitted
    text(max(x)/2,max(y)*0.8,['\alpha=' num2str(-slope,'%0.1f')],'color',ColorDark,'FontSize',14,'FontWeight','bold')
    text(max(x)/2,max(y)*0.6,['R^2=' num2str(R2,'%0.3f')],'color',ColorDark,'FontSize',14,'FontWeight','bold')
    
    % 6. Plot loglog peaks histogramm
    axes('outerposition',aHistPeaksLog)
    loglog(x,y,'.','color',ColorLight,'MarkerSize',30)
    hold on
    loglog(xfit,yfit,'--','color',ColorDark,'linewidth',2.5)
    set(gca,'box','off','FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    xlabel('# peaks')
    
    
    % Plot peaks groups in raster
    Data=Experiment.Data.Data;
    F=Experiment.Data.Frames;
    C=Experiment.Data.Neurons;
    Ac=Experiment.Raster.Activity;
    Co=Experiment.Raster.Coactivity;
    ThSync=Experiment.Peaks.Threshold;
    
    peaks=max(PeaksIdx);
    fpm=str2num(get(handles.FPSEdit,'string'))*60;
    
    Experiment=PlotSimpleExperiment(Experiment, handles);
    % Raster
    idx_vector_state=Experiment.Clustering.VectorStateIndex;
    axes(findobj('tag',['RasterAxes' DataName])); hold on
    for i=1:peaks
        peak=find(PeaksIdx==i);
        group=idx_vector_state(i);
        for j=1:length(peak)
            idxActive=find(Data(peak(j),:));
            plot(peak(j),idxActive,'.','color',GroupColors(group,:),'MarkerSize',15)
        end
    end
    
    % Plot peaks
    axes(findobj('tag',['CoactiveAxes' DataName])); hold on
    plot(get(gca,'xlim'),[ThSync-.5 ThSync-.5],'--','color',[0 0 0],'lineWidth',2)
    for i=1:peaks
        peak=find(PeaksIdx==i);
        group=idx_vector_state(i);
        for j=1:length(peak)
            x=[peak(j)-5:peak(j)+5];
            idx=find(x<=0);
            x(idx)=1;
            idx=find(x>=F);
            x(idx)=F;
            plot(x/fpm,Co(x),'-','color',GroupColors(group,:),'linewidth',1.5)
        end
    end
    
    % Plot Activity
    axes(findobj('tag',['ActivityAxes' DataName])); hold on
    for i=1:C
        ActiveCell=find(Active==i);
        if ActiveCell
            color=ColorComb(IndexesCombClean(ActiveCell),:);
        else
            color=[0 0 0];
        end
        plot(i,Ac(i),'.','color',color,'MarkerSize',15,'linewidth',2)
    end
    
    
    % Plot cells gruop in raster
    C=size(Data,2);
    AllCellsGroups=zeros(Groups,C);
    AllCellsGroups(:,Active)=CellsGroups;
    figure(findobj('name',['Experiment (' DataName ')']))
    %figure('name','Groups','position',[0 0 1220 460])
    axes('outerposition',[0 .34 .125 .66]);
    image(AllCellsGroups'+1)
    colormap([1 1 1; GroupColors])
    set(gca,'Tag',['GroupAxesRaster' DataName])
    set(gca,'YDir','normal','xcolor',ColorTextLines,'ycolor',ColorTextLines,'XTick',1:Groups)
    set(gca,'Ytick',[])
    xlabel('State'); ylabel('# neuron')
    title('States')
    hold on
    for i=1:Groups-1
        plot([i+.5 i+.5],[0 C+2],'-k')
    end
    for i=1:C-1
        plot([0 Groups+0.5],[i+0.5 i+0.5],'-k')
    end
    
    % Outputs
    Experiment.Ranks.States=IndexesClean;
    Experiment.Network.PeaksSlope=-slope;
end

%% Plot Adjacency Matrix
function Experiment=PlotAdjacent(Experiment,handles)
    
    Data=Experiment.Data.Data;
    DataName=Experiment.Data.DataName;
    Active=Experiment.Peaks.ActiveNeurons;
    ActiveDataPeaks=Experiment.Peaks.ActiveDataPeaks;
    
    obj=findobj('name',['Correlation (' DataName ')']);
    if isempty(obj)
        figure('name',['Correlation (' DataName ')'],'numbertitle','off','position',[0 0 900 600]);
    else
        figure(obj);
    end
    clf
    
    % Colors
    % BGColor=[0.7 0.7 0.8];
    BGColor=[1 1 1];
    ColorLight=[0.9 0 0.7];
    ColorDark=[0.5 0 0.3];
    ColorTextLines=[0 0 0];
    
    % Axes
    aCorr=[0 0.5 0.33 0.5];
    aAdjacent=[.33 0.5 0.33 .5];
    aLinks=[0.66 0.5 0.34 0.5];
    aHistLinks=[0 0 0.5 0.5];
    aHistLinksPower=[0.5 0 0.5 0.5];
    
    % others
    aCent=[0.25 0.33 0.375 0.33];
    aHistCent=[0.625 0.33 0.375 0.33];
    aGeo=[0.25 0 0.375 0.33];
    aHistGeo=[0.625 0 0.375 0.33];
    
    % Plot correlation
    axes('outerposition',aCorr)
    set(gca,'Tag',['ImgCorrelationAxes' DataName])
    
    % Get Correlation method
    %{
    Num=get(handles.CorrMethodPopupmenu,'Value');
    Strings=get(handles.CorrMethodPopupmenu,'String');
    Method=Strings{Num};
    if strcmp(Method,'Jaccard')
        ImgCorr=1-squareform(pdist(Data',Method));
    else
        ImgCorr=corr(Data,'type',Method);
    end
    %}
    %***************
    % Adjacency matrix from peaks
    C=size(Data,2);
    ImgCorr=GetAdjacencyFromPeaks(handles,ActiveDataPeaks); % just neurons in peaks
    %ImgCorr=(Data()'*Data()).*(1-eye(C)); % All neurons
    
    % JP modified oct-13
    AdjW=zeros(C);
    AdjW(Active,Active)=ImgCorr;
    imagesc(AdjW)

    set(gca,'YDir','normal','Xtick',[10:10:C],...
        'Ytick',[10:10:C],'xcolor',ColorTextLines,'ycolor',ColorTextLines)
    %title([Method ' correlation'],'color',ColorTextLines)
    title(['Adjacency matrix (weighted)'],'color',ColorTextLines)
    xlabel('# neuron')
    
       
    axes('outerposition',aAdjacent)
    ThCorr=get(handles.CorrThSlider,'value')*100;
    AdjB=AdjW>=ThCorr; % >= Threshold because of WeightTh is more than 
    imagesc(AdjB)
    
    set(gca,'YDir','normal','Xtick',[10:10:C],...
        'Ytick',[10:10:C],'xcolor',ColorTextLines,'ycolor',ColorTextLines)
    set(gca,'Tag',['ImgAdjacentCorrAxes' DataName])
    title([' Adjacency matrix'],'color',ColorTextLines)
    xlabel('# neuron')
    
    % Plot links
    axes('outerposition',aLinks)
    Links=sum(AdjB,1);
    
    
    [value Experiment.Ranks.Links]=sort(Links,'descend');
    if ~Links
        Experiment.Adjacency.AdjacentMatrix=ImgCorr;
        Experiment.Adjacency.AdjacentMatrixThreshold=AdjB;
        Experiment.Adjacency.Links=Links;
        return;
    end
    %plot(Active,Links,'-k')
    plot(Links,'-k')
    set(gca,'xlim',[0 C],'ylim',[0 max(Links)+1])
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    title('Links','color',ColorTextLines)
    xlabel('# neuron')
    ylabel('Count of links')
    set(gca,'Tag',['LinksAxes' DataName])

    % Plot links histogram
    idx=Links>0;
    Links=Links(idx);
    axes('outerposition',aHistLinks)
    max_Links=max(Links);
    if max_Links
        N=length(Links);
        nbins=round(log2(N)+1);
        [y x]=hist(Links,nbins);%This is a probe ,0:max_Links);
        idx=find(y);
        y=y(idx);
        x=x(idx);
        
        
        plot(x,y,'.','color',ColorLight,'MarkerSize',30)
        set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
        %[y x]=hist(Links,0:max_Links);
        %bar(x(2:max_Links+1),y(2:max_Links+1),'facecolor',ColorLight,'edgecolor',ColorDark)
        hold on
        set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
        %title('Links histogram','color',ColorTextLines)
        xlabel('Links')
        ylabel('Count of neurons')
    end
    set(gca,'Tag',['LinksHistAxes' DataName])

    % Plot power law fit (only neurons with at least one link)
    [slope, intercept, R2] = logfit(x,y,'loglog');
    xfit=min(x):(max(x)-min(x))/100:max(x);
    yfit=(10^intercept)*xfit.^(slope);
    plot(xfit,yfit,'--','color',ColorDark,'linewidth',2.5)
    % Write the function fitted
    text(max(x)/2,max(y)*0.8,['\alpha=' num2str(-slope,'%0.1f')],'color',ColorDark,'FontSize',14,'FontWeight','bold')
    text(max(x)/2,max(y)*0.6,['R^2=' num2str(R2,'%0.3f')],'color',ColorDark,'FontSize',14,'FontWeight','bold')

    % Loglog activity histogram
    axes('outerposition',aHistLinksPower)
    loglog(x,y,'.','color',ColorLight,'MarkerSize',30)
    hold on
    loglog(xfit,yfit,'--','color',ColorDark,'linewidth',2.5)
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    set(gca,'box','off')
     
    % Network properties
    %{
    obj=findobj('name',['Network properties (' DataName ')']);
    if isempty(obj)
        figure('name',['Network properties (' DataName ')'],'numbertitle','off','position',[0 0 900 300]);
    else
        figure(obj);
    end
    clf
    
    % Colors
    BGColor=[1 1 1];
    ColorLight=[0.9 0 0.7];
    ColorDark=[0.5 0 0.3];
    ColorTextLines=[0 0 0];
    
    % Axes
    aCent=[0 0 0.5 1];
    aGeo=[0.5 0 0.5 1];
    
    % Plot centrality (changes in minimum path)
    axes('outerposition',aCent)
    %}
    
    %[ChangesInMinPaths TrimsInNet]=Centrality_JP(AdjB(Active,Active));
    %plot(Active,ChangesInMinPaths,'-k')
    %[value Experiment.Ranks.Centrality]=sort(ChangesInMinPaths,'descend');
    
    %{
    
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    title('Centrality','color',ColorTextLines)
    xlabel('# neuron')
    ylabel('Count of changes')
    set(gca,'Tag',['ChangesAxes' DataName])
    %legend('Changes in minimum paths','Trims in network')
    
    % Plot Geodesic
    %
    axes('outerposition',aGeo)
    [Geo GeoSingleMean GeoMean Trims]=Geodesic_JP(AdjacentCorr);
    plot([0 C],[GeoMean GeoMean],'--k')
    hold on
    plot(Active,mean(Geo),'ok')
    title('Neuron geodesic distance');xlabel('# neuron'); ylabel('Geodesic mean')
    set(gca,'box','off')
    %}
    
    % Get cardinals
    Ca=length(AdjB(Active,Active));
    Cardinality=zeros(Ca,1);
    for i=1:Ca
        Cardinality(i)=find(Experiment.Ranks.Activity==i);
        Cardinality(i)=Cardinality(i)+find(Experiment.Ranks.Peaks==i);
        Cardinality(i)=Cardinality(i)+find(Experiment.Ranks.Links==i);
%        Cardinality(i)=Cardinality(i)+find(Experiment.Ranks.Centrality==i);
    end
    [value Experiment.Ranks.Cardinality]=sort(Cardinality,'ascend');

    % Outputs
    Experiment.Adjacency.AdjacentMatrix=ImgCorr;
    Experiment.Adjacency.AdjacentMatrixThreshold=AdjB;
    Experiment.Adjacency.Links=sum(AdjB,1);
%    Experiment.Adjacency.Centrality=ChangesInMinPaths;
    Experiment.Network.LinksSlope=-slope;
    
end

%% Plot Network
function Experiment=PlotNetwork(Experiment,handles)

    % Inputs
    DataName=Experiment.Data.DataName;
    Active=Experiment.Peaks.ActiveNeurons;
    Adjacency=Experiment.Adjacency.AdjacentMatrix;
    Groups=Experiment.Clustering.TotalStates;
    IndexesComb=Experiment.Structure.StatesCombinations;
    IdxGroups=Experiment.Ranks.States;
    GroupColors=Experiment.Colors.States;
    C=Experiment.Data.Neurons;
    
    % Figure Network
    obj=findobj('name',['Network (' DataName ')']);
    if isempty(obj)
        figure('name',['Network (' DataName ')'],'numbertitle','off','position',[0 0 600 600]);
    else
        figure(obj);
    end
    clf
    axes('outerposition',[0 0 1 1])
    hold on
    
    % Colors
    BGColor=[1 1 1];
        
    % Combinations
    Comb=2^Groups-1;
    ColorComb=zeros(Comb,3);
    Level{Groups}=[];
    for i=1:Comb
        bin=dec2bin(i,Groups);
        color=[];
        for j=1:Groups
            if str2num(bin(j))
                color=[color; GroupColors(j,:)];
            end
        end
        ColorComb(i,:)=mean(color,1);

        % Levels
        for j=1:Groups
            if size(color,1)==j
                if isempty(Level{j})
                    Level{j}=find(IndexesComb==i);
                else
                    Level{j}=[Level{j} 0 find(IndexesComb==i)];%0 identifies different groups inside of a level
                end
            end
        end
    end
    
    XY=Experiment.Data.XY;
    XYsort=get(handles.XYSortPopupmenu,'value');
    
    C=length(Adjacency);
    if isempty(XY)
        % XY if there are not coordinates
        XY_circle=1;
        switch (XYsort)
            case 1 % None
                XY=[cos(2*pi*[1:C]'/C) sin(2*pi*[1:C]'/C)];
            case 3 % Activity
                idx=Experiment.Ranks.Activity;
                XY(idx,:)=[cos(2*pi*[1:C]'/C) sin(2*pi*[1:C]'/C)];
            case 4 % Peaks
                idx=Experiment.Ranks.Peaks;
                XY(idx,:)=[cos(2*pi*[1:C]'/C) sin(2*pi*[1:C]'/C)];
            case 5 % Links
                idx=Experiment.Ranks.Links;
                XY(idx,:)=[cos(2*pi*[1:C]'/C) sin(2*pi*[1:C]'/C)];
            case 6 % Centrality
                idx=Experiment.Ranks.Centrality;
                XY(idx,:)=[cos(2*pi*[1:C]'/C) sin(2*pi*[1:C]'/C)];
            case 7 % Cardinality
                idx=Experiment.Ranks.Cardinality;
                XY(idx,:)=[cos(2*pi*[1:C]'/C) sin(2*pi*[1:C]'/C)];
            case 2 % Groups
                g2=0;
                for i=1:Groups
                    g2=g2+length(find(Level{:,i}==0));
                end
                CsLs=Groups+g2;
                C2=C+CsLs*5; % Original
                C2=C+7*5; % Ad hoc
                XY=zeros(C,2);
                Ci=1;
                for i=1:Groups
                    CL=length(Level{i});
                    if CL
                        subLevel=[0 Level{i} 0];
                        idxk=find(subLevel==0);
                        k=length(idxk)-1;
                        for j=1:k
                            CsL=idxk(j+1)-idxk(j)-1;
                            XYL=[cos(2*pi*[Ci:(Ci+CsL-1)]'/C2) sin(2*pi*[Ci:(Ci+CsL-1)]'/C2)];
                            idxsL=subLevel((idxk(j)+1):(idxk(j+1)-1));
                            XY(idxsL,:)=XYL;
                            Ci=Ci+CsL+5; % Original
                        end
                    end
                end
        end
        set(gca,'xlim',[-1.2 1.2],'ylim',[-1.2 1.2])
    else
        set(handles.XYSortPopupmenu,'value',1);
        XY_circle=0;
        XY=XY(Active,:);
        
        % Normalize
        %{
        Xmax=max(XY(:,1));
        Ymax=max(XY(:,2));
        Xmin=min(XY(:,1));
        Ymin=min(XY(:,2));
        XY(:,1)=XY(:,1)-Xmin;
        XY(:,2)=XY(:,2)-Ymin;
        dx=Xmax-Xmin;
        dy=Ymax-Ymin;
        if dx>dy
            dmax=dx;
        else
            dmax=dy;
        end
        XY(:,1)=XY(:,1)/dmax;
        XY(:,2)=XY(:,2)/dmax;
        set(gca,'xlim',[-0.1 1.1],'ylim',[-0.1 1.1])
        %}
        set(gca,'xlim',[0 512],'ylim',[0 512],'xcolor',[1 1 1],'ycolor',[1 1 1])
    end
    Experiment.Network.XY=XY;
    
    % Plot edges
    y=squareform(Adjacency,'tovector');
    ThCorr=get(handles.CorrThSlider,'value')*100;
    idx=y>ThCorr;
    y=sort(y(idx));
    Q1=median(y(find(y<median(y))));
    Q2=median(y); 
    Q3=median(y(find(y>median(y))));
    XYedgecolor=[];
    %width=2;
    %
    for a=1:C-1
        for b=a+1:C
            if Adjacency(a,b)>=ThCorr
                color=mean([ColorComb(IndexesComb(a),:);...
                    ColorComb(IndexesComb(b),:)]);
                %color=[0 0 0];
                if Adjacency(a,b)< Q1
                    width=1;
                else
                    width=2;
                end
                plot(XY([a b],1),XY([a b],2),'color',color,'LineWidth',width);
                XYedgecolor(a,b,:)=color;
                XYedgecolor(b,a,:)=color;
                hold on
            end
        end
    end
    %}
    Experiment.Network.XYedgecolor=XYedgecolor;
    
    % Plot neurons
    for i=1:C
        color=ColorComb(IndexesComb(i),:);
        %color=[1 1 1];
        plot(XY(i,1),XY(i,2),'.k','MarkerSize',35)
        plot(XY(i,1),XY(i,2),'.','color',color,'MarkerSize',23)
        %plot(XY(i,1),XY(i,2),'o','color',color,'linewidth',2,'markersize',20)
        XYnodecolor(i,:)=color;
        if XY_circle
            text(XY(i,1)*1.1,XY(i,2)*1.1,num2str(Active(i)),'HorizontalAlignment','Center')
        else
            text(XY(i,1)+0.02,XY(i,2)-0.02,num2str(Active(i)),'HorizontalAlignment','Center')
        end
    end
    set(gca,'view',[0 -90])
    set(gca,'xtick',[],'ytick',[])
    title(['thin lines < ' num2str(Q1)])
    
    Experiment.Network.XYnodecolor=XYnodecolor;
    
    
    % Network properties
    Adjacency=Adjacency>=ThCorr;
    [Adjacency idxConnected idxNotConnected] = OnlyConnected(Adjacency);
    
    Experiment.Network.idxConnected=idxConnected;
    Experiment.Network.idxNotConnected=idxNotConnected;
    if length(idxNotConnected)
        % idx
        IdxCard=Experiment.Ranks.Cardinality;
        j=1;
        for i=idxNotConnected
            id(j)=find(IdxCard==i);
            j=j+1;
        end
        idxActive=setdiff(1:(max(IdxCard)),id);
        IdxCard=IdxCard(idxActive);

        [temp idxSort]=sort(IdxCard);
        IdxCard(idxSort)=1:length(idxConnected);

        Experiment.Ranks.Cardinality=IdxCard;
    end
    Experiment.Adjacency.AdjacencyOnlyConnected=Adjacency;
    Plot_Network_Properties(Experiment);
    
    % Plot Neurons features
    obj=findobj('name',['Neurons features (' DataName ')']);
    if isempty(obj)
        figure('name',['Neurons features (' DataName ')'],'numbertitle','off','position',[0 0 600 600]);
    else
        figure(obj);
    end
    clf
    
    % Colors
    BGColor=[1 1 1];
    ColorLight=[0.9 0 0.7];
    ColorDark=[0.5 0 0.3];
    ColorTextLines=[0 0 0];
    
    % Axes
    aPeaks=[0 0.75 1 0.25];
    aActivity=[0 0.5 1 0.25];
    aLinks=[0 0.25 1 0.25];
    aCent=[0 0 1 0.25];
    
    % Get data
    CellsActGroups=Experiment.Structure.RateNeuronsStates;
    Peaks=max(Experiment.Peaks.Index);
    IndexesCombClean=Experiment.Structure.StatesCombinations;
    maxCellPeaks=max(sum(CellsActGroups,1))/Peaks*100;
    Activity=Experiment.Raster.Activity;
    Links=Experiment.Adjacency.Links;
%    ChangesInMinPaths=Experiment.Adjacency.Centrality;
    C=length(Active);
    N=Experiment.Data.Neurons;
    
    % 1. Plot peaks
    axes('outerposition',aPeaks); hold on
    for i=1:C
        bar(Active(i),sum(CellsActGroups(:,i))/Peaks*100,'facecolor',ColorComb(IndexesCombClean(i),:),'edgecolor',[0 0 0])
    end
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    set(gca,'Xlim',[.5 N+.5],'XTick',[])
    set(gca,'Ylim',[0 maxCellPeaks])
    ylabel('% peaks'); title('Peaks')
    set(gca,'Tag',['PeaksFAxes' DataName])
    
    % 2. Plot activity
    axes('outerposition',aActivity); hold on
    bar([1:length(Activity)],Activity,'facecolor',[1 1 1],'edgecolor',[0 0 0])
    for i=1:C
        bar(Active(i),Activity(Active(i)),'facecolor',ColorComb(IndexesCombClean(i),:),'edgecolor',[0 0 0])
    end
    set(gca,'Xlim',[.5 N+.5],'xtick',[],'ylim',[0 max(Activity)])
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    title('Activity','color',ColorTextLines)
    ylabel('% active frames')
    
    % 3. Plot links
    axes('outerposition',aLinks); hold on
    for i=1:C
        bar(Active(i),Links(Active(i)),'facecolor',XYnodecolor(i,:),'edgecolor',[0 0 0])
    end
    if (max(Links))
        set(gca,'Xlim',[.5 N+.5],'xtick',[],'ylim',[0 max(Links)])
        set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
        title('Links','color',ColorTextLines)
        set(gca,'Tag',['LinksFAxes' DataName])
        ylabel('Count of links')
    else
        return;
    end
    
    % 4. Plot local clustering coefficient
    Clocal=clustering_coef_bu(Adjacency);
    axes('outerposition',aCent); hold on
    C=length(Clocal);
    for i=1:C
        bar(Active(idxConnected(i)),Clocal(i),'facecolor',XYnodecolor(idxConnected(i),:),'edgecolor',[0 0 0])
    end
    set(gca,'Xlim',[.5 N+.5],'ylim',[0 1])
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    title('Local clustering coefficient','color',ColorTextLines)
    ylabel('C_{local}')
    xlabel('# neuron')
    set(gca,'Tag',['LCCFAxes' DataName])
    
    % Plot Hierarchy
    %Plot_Hierarchy(Experiment);
    Experiment=Plot_Hierarchy_Color(Experiment);
    %Experiment=Plot_Hierarchy_Mean(Experiment);
    Experiment=Plot_DegreeDistribution_Color(Experiment);
    Experiment=SmallWorld_JP3(Experiment);
    

    % Neurons to begin disconnection
    %Experiment.Network.Disconnection=find(Cuts2,1,'first')/N;
    
    % Outputs
    %{
    Experiment.Network.Assortativity=A;
    Experiment.Network.GeodesicMean=Gmean;
    Experiment.Network.Transitivity=T;
    Experiment.Network.Connection=LinksRatio;
    %}
end

%% ...::: Create Objects :::...

function GetFilePushbutton_CreateFcn(hObject, eventdata, handles)
end
function AnalizePushbutton_CreateFcn(hObject, eventdata, handles)
end
function ProcessText_CreateFcn(hObject, eventdata, handles)
    set(hObject,'String',{''})
end
%% Data Panel
function DataUipanel_CreateFcn(hObject, eventdata, handles)
end
function PlotDataPushbutton_CreateFcn(hObject, eventdata, handles)
end
function DataText_CreateFcn(hObject, eventdata, handles)
end
function DataPopupmenu_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',{'---None---'})
end
function XYText_CreateFcn(hObject, eventdata, handles)
end
function XYPopupmenu_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',{'---None---'})
end
function FPSEdit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String','4')
end

%% Activity Panel
function ActivityUipanel_CreateFcn(hObject, eventdata, handles)
end
function GetActPushbutton_CreateFcn(hObject, eventdata, handles)
end

%% Peaks Panel
function PeaksUipanel_CreateFcn(hObject, eventdata, handles)
end
function PlotPeaksPushbutton_CreateFcn(hObject, eventdata, handles)
end
function PeaksThText_CreateFcn(hObject, eventdata, handles)
end
function PeaksThSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    set(hObject,'Value',6)
end
function PeaksThEdit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',6)
end
function MCPushbutton_CreateFcn(hObject, eventdata, handles)
end
function JoinPeaksCheckbox_CreateFcn(hObject, eventdata, handles)
    set(hObject,'Value',0)
end
function BinaryPeaksCheckbox_CreateFcn(hObject, eventdata, handles)
    set(hObject,'Value',0)
end
function SimMethodText_CreateFcn(hObject, eventdata, handles)
end
function SimMethodPopupmenu_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    DataStrings={'Cosine','Euclidean','Jaccard','LLEdistance','Seuclidean','Mahalanobis',...
        'Cityblock','Minkowski','Correlation','Spearman','Hamming','Chebychev'};
    set(hObject,'String',DataStrings)
end

%% Dimensional Reduction Panel
function DimRedUipanel_CreateFcn(hObject, eventdata, handles)
end
function NoRedRadiobutton_CreateFcn(hObject, eventdata, handles)
end
function RedRadiobutton_CreateFcn(hObject, eventdata, handles)
end
function DimRedText_CreateFcn(hObject, eventdata, handles)
end
function DimRedPopupmenu_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    DataStrings={'LLE','PCA'};
    set(hObject,'String',DataStrings)
end
function NeighborsText_CreateFcn(hObject, eventdata, handles)
end
function NeighborsEdit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',1)
end

%% Clustering Panel
function HierarchicalRadiobutton_CreateFcn(hObject, eventdata, handles)
    set(hObject,'value',1)
end
function ClustPopupmenu_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    DataStrings={'Average','Centroid','Complete','Median',...
        'Single','Ward','Weighted'};
    set(hObject,'String',DataStrings)
end
function GroupsThSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    set(hObject,'String',3)
end
function GroupsThEdit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',3)
end
function CleanThSlider_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
    set(hObject,'Value',75)
end
function CleanThEdit_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',75)
end

%% Network Panel
function CorrMethodPopupmenu_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    DataStrings={'Synchronization','Jaccard','Pearson','Kendall','Spearman'};
    set(hObject,'String',DataStrings)
end
function XYSortPopupmenu_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    DataStrings={'---None---','Groups','Activity','Peaks','Links','Centrality','Cardinality'};
    set(hObject,'String',DataStrings)
end

%% Colors Panel
function GroupColorsPopupmenu_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',{'R G B Y M C O P Ro Br'})
end
function XYColorsPopupmenu_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
    set(hObject,'String',{'Black'})
end

%% ...::: Callbacks :::...

%% Menu
function DataMenu_Callback(hObject, eventdata, handles)
end
function SelectFileMenu_Callback(hObject, eventdata, handles)
    [FileName,PathName] = uigetfile('*.mat','Select the EXPERIMENT file');
    if isequal(FileName,0)
        %disp('User selected Cancel')
    else
        %evalin('base', 'clear;')
        temp_data=load(fullfile(PathName,FileName)); % will be loaded as structure
        file_variables=fieldnames(temp_data);% get the field names as cell array
        for ii=1:length(file_variables)
            % file_variables{ii} - string of the field name
            % temp_data.(file_variables{ii}) - dynamic field reference
        assignin('base', file_variables{ii}, temp_data.(file_variables{ii}));
        end       
        %Processing('start',handles)
            DataStrings = evalin('base','who');
            set(handles.DataPopupmenu,'String',[{'---None---'};DataStrings])
            set(handles.DataPopupmenu,'Value',1)
            set(handles.XYPopupmenu,'String',[{'---None---'};DataStrings])
            set(handles.XYPopupmenu,'Value',1)
            set(handles.GroupColorsPopupmenu,'String',[{'R G B Y M C O P Ro Br'};DataStrings])
            set(handles.GroupColorsPopupmenu,'Value',1)
            set(handles.XYColorsPopupmenu,'String',[{'Black'};DataStrings])
            set(handles.XYColorsPopupmenu,'Value',1)
        %Processing('end',handles)
    end   
end
function AnalizeMenu_Callback(hObject, eventdata, handles)
    %Processing('start',handles)
    Update(handles);
    %Processing('end',handles)
end
function ExitMenu_Callback(hObject, eventdata, handles)
    close(gcbf);
end
function HelpMenu_Callback(hObject, eventdata, handles)
end
function AboutMenu_Callback(hObject, eventdata, handles)
    msg=sprintf(['Toolbox for neural network analysis.'...
        '\nBy Pérez-Ortega Jesús Esteban.\nVersion 5.6\nAug - 2014']);
    msgbox(msg,'About')
end
function AnalizePushbutton_Callback(hObject, eventdata, handles)
    %Processing('start',handles)
    Update(handles);
    %Processing('end',handles)
end
function LoadDataPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadDataPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    DataStrings = evalin('base','who');
    set(handles.DataPopupmenu,'String',[{'---None---'};DataStrings])
    set(handles.DataPopupmenu,'Value',1)
    set(handles.XYPopupmenu,'String',[{'---None---'};DataStrings])
    set(handles.XYPopupmenu,'Value',1)
    set(handles.GroupColorsPopupmenu,'String',[{'R G B Y M C O P Ro Br'};DataStrings])
    set(handles.GroupColorsPopupmenu,'Value',1)
    set(handles.XYColorsPopupmenu,'String',[{'Black'};DataStrings])
    set(handles.XYColorsPopupmenu,'Value',1)
end

%% Data Panel
function PlotDataPushbutton_Callback(hObject, eventdata, handles)
    global green_image step_execution gray_image waiting_image
    %Processing('start',handles)
    Processing('msj','Plotting Simple Experiment',handles)
    Experiment=Get_Experiment(handles);
    if isfield(Experiment,'Data')
        if isfield(Experiment.Data,'Data')     
            set(handles.PlotDataPushbutton,'CData', waiting_image);
            set(handles.PlotPeaksPushbutton,'CData', gray_image);
            set(handles.DimRedPushbutton,'CData', gray_image);
            set(handles.ClusteringPushbutton,'CData', gray_image);
            set(handles.PlotNetworkPushbutton,'CData', gray_image);
            Experiment=PlotSimpleExperiment(Experiment,handles);
            assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
            Experiment=PlotActivity(Experiment,handles);
            assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
            set(handles.PlotDataPushbutton,'CData', green_image);
            step_execution = [1 0 0 0 0];
            Processing('end','Plotting Simple Experiment',handles)
        else
            Processing('error','Plotting Simple Experiment',handles)
        end
    else
        Processing('error','Getting Simple Experiment',handles)
    end
    %Processing('end',handles)
end
function DataPopupmenu_Callback(hObject, eventdata, handles)
    global gray_image step_execution
    % Get data
    set(handles.PlotDataPushbutton,'CData', gray_image);
    set(handles.PlotPeaksPushbutton,'CData', gray_image);
    set(handles.DimRedPushbutton,'CData', gray_image);
    set(handles.ClusteringPushbutton,'CData', gray_image);
    set(handles.PlotNetworkPushbutton,'CData', gray_image);
    step_execution = [0 0 0 0 0];
    
    Processing('msj','Reading Experiment',handles)
    Experiment=Get_Data(handles);
    if isempty(Experiment.Data.Data)
        Processing('error','Reading Experiment',handles)
        return;
    end
    Processing('msj','Getting AcCoMax',handles)
    [AcMax CoMax]=GetAcCoMax(Experiment.Data.Data);
    Processing('msj','Getting CoTh',handles)
    CoTh=get(handles.PeaksThSlider,'Value');
    if CoTh>CoMax
        CoTh=CoMax;
    end
    set(handles.PeaksThSlider,'Value',CoTh);
    set(handles.PeaksThSlider,'Max',CoMax)
    set(handles.PeaksThSlider,'SliderStep',1./[CoMax CoMax/2])
    set(handles.PeaksThEdit,'String',num2str(CoTh))
    
    % Get other data
    Processing('msj','Getting Other Data',handles)
    Experiment=Get_XY(Experiment,handles);
    Experiment=Get_XYColors(Experiment,handles);
    Experiment=Get_GroupColors(Experiment,handles);
    
    % Set to workspace
    Processing('msj','Setting to Workspace',handles)
    assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
    Processing('end','Selected Experiment',handles)
    %Processing('end',handles)
end

function XYPopupmenu_Callback(hObject, eventdata, handles)
    DataPopupmenu_Callback(hObject, eventdata, handles)
end
function FPSEdit_Callback(hObject, eventdata, handles)
    Processing('msj','Getting Frames per Second',handles)
    user_entry = str2double(get(hObject,'string'));
    if isnan(user_entry)
        set(hObject,'string',4);
    elseif(user_entry>10)
        set(hObject,'string',10);
    elseif(user_entry<=0)
        set(hObject,'string',0.1);
    else
        set(hObject,'string',num2str(user_entry,'%.1f'));
    end
    DataPopupmenu_Callback(hObject, eventdata, handles)
end

%% Activity Panel
function GetActPushbutton_Callback(hObject, eventdata, handles)
    global gray_image step_execution
    set(handles.PlotDataPushbutton,'CData', gray_image);
    set(handles.PlotPeaksPushbutton,'CData', gray_image);
    set(handles.DimRedPushbutton,'CData', gray_image);
    set(handles.ClusteringPushbutton,'CData', gray_image);
    set(handles.PlotNetworkPushbutton,'CData', gray_image);
    step_execution = [0 0 0 0 0];
    Processing('msj','Getting Activity Distribution',handles)
    Experiment=Get_Experiment(handles);
    if isfield(Experiment,'Raster')
        Processing('msj','Plotting Activity Distribution',handles)
        Experiment=PlotActivity(Experiment,handles);
        assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
        Processing('end','Plotting Activity Distribution',handles)
    else
        Processing('error','Plotting Activity Distribution',handles)
    end
    %Processing('end',handles)
    
end

%% Peaks Panel
function PlotPeaksPushbutton_Callback(hObject, eventdata, handles)
    global green_image step_execution yellow_image gray_image waiting_image
    if step_execution(1)==1
        set(handles.PlotPeaksPushbutton,'CData', waiting_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        Processing('msj','Plotting Peaks',handles)
        Experiment=Get_Experiment(handles);
        if isfield(Experiment,'Raster')
            Experiment=PlotPeaks(Experiment,handles);
            assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
            set(handles.PlotPeaksPushbutton,'CData', green_image);
            
            step_execution = [1 1 0 0 0];
            Processing('end','Plotting Peaks',handles)
        else
            Processing('error','Plotting Peaks',handles)
        end
        %Processing('end',handles)
   else
        set(handles.PlotDataPushbutton,'CData', yellow_image);
        set(handles.PlotPeaksPushbutton,'CData', yellow_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [2 2 0 0 0];
        Processing('error','Incorrect Step',handles)
    end
end
function PeaksThSlider_Callback(hObject, eventdata, handles)
    global gray_image step_execution green_image 
    if step_execution(1)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 0 0 0 0];
    else
        set(handles.PlotDataPushbutton,'CData', gray_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [0 0 0 0 0];
    end
    
    slider_value = round(get(handles.PeaksThSlider,'Value'));
    set(handles.PeaksThEdit,'String',slider_value);
end
function PeaksThEdit_Callback(hObject, eventdata, handles)
global gray_image step_execution green_image
    if step_execution(1)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 0 0 0 0];
    else
        set(handles.PlotDataPushbutton,'CData', gray_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [0 0 0 0 0];
    end
    user_entry = str2double(get(hObject,'string'));
    if isnan(user_entry)
        PeaksTh=get(handles.PeaksThSlider,'Value');
        set(hObject,'string',num2str(PeaksTh));
        return;
    end
    set(hObject,'string',floor(user_entry));
        
    PeaksMax=get(handles.PeaksThSlider,'Max');
    if(user_entry>PeaksMax)
        set(hObject,'string',num2str(PeaksMax));
        set(handles.PeaksThSlider,'Value',PeaksMax);
    elseif(user_entry<0)
        set(hObject,'string','0');
        set(handles.PeaksThSlider,'Value',0);
    else
        set(handles.PeaksThSlider,'Value',user_entry);
    end
end
function MCPushbutton_Callback(hObject, eventdata, handles)
    global gray_image step_execution green_image
    if step_execution(1)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 0 0 0 0];
    else
        set(handles.PlotDataPushbutton,'CData', gray_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [0 0 0 0 0];
    end
    Experiment=Get_Experiment(handles);
    if isfield(Experiment,'Raster')
        options.Resize='on';
        answer = inputdlg('Enter number of iterations:','Montecarlo test',...
        1,{'1000'},options);
        if isempty(answer)
            %Processing('end',handles)
            return;
        else
            iterations=str2num(answer{1});
            if ~isempty(iterations)
                if iterations<2
                    iterations=2;
                end
                Data=Experiment.Data.Data;
                DataName=Experiment.Data.DataName;
                obj=findobj('name',['Montecarlo Test (' DataName ')']);
                if isempty(obj)
                    obj=figure('name',['Montecarlo Test (' DataName ')'],'numbertitle','off','position',[0 0 600 300]);
                else
                    figure(obj);
                end
                clf
                tic
                [th P]=MC_JP(Data,iterations,0.05,obj);
                toc
                Experiment.Test.Montecarlo=P;
                assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
            end
        end
    end    
    %Processing('end',handles)
end
function JoinPeaksCheckbox_Callback(hObject, eventdata, handles)
    global gray_image step_execution green_image
    if step_execution(1)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 0 0 0 0];
    else
        set(handles.PlotDataPushbutton,'CData', gray_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [0 0 0 0 0];
    end
end
function BinaryPeaksCheckbox_Callback(hObject, eventdata, handles)
    global gray_image step_execution green_image
    if step_execution(1)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 0 0 0 0];
    else
        set(handles.PlotDataPushbutton,'CData', gray_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [0 0 0 0 0];
    end
end
function SimMethodPopupmenu_Callback(hObject, eventdata, handles)
    global gray_image step_execution green_image
    if step_execution(1)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 0 0 0 0];
    else
        set(handles.PlotDataPushbutton,'CData', gray_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [0 0 0 0 0];
    end
end

%% Dimensional Reduction Panel
function DimRedPopupmenu_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 0 0 0];
    end
end
function DimRedPushbutton_Callback(hObject, eventdata, handles)
    global green_image step_execution yellow_image gray_image waiting_image
    if step_execution(1)==1 && step_execution(2)==1
        set(handles.DimRedPushbutton,'CData', waiting_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        Processing('msj','Dimensional Reduction',handles)
        Experiment=Get_Experiment(handles);
        if isfield(Experiment,'Peaks')
            Experiment=PlotDimRed(Experiment,handles);
            assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
            set(handles.DimRedPushbutton,'CData', green_image);            
            step_execution = [1 1 1 0 0];
            Processing('end','Dimensional Reduction',handles)
        else
            Processing('error','Dimensional Reduction',handles)
        end
        %Processing('end',handles)
    elseif step_execution(1)==1 && (step_execution(2)==0 || step_execution(2)==2)
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', yellow_image);
        set(handles.DimRedPushbutton,'CData', yellow_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 2 2 0 0];
        Processing('error','Incorrect Step',handles)
   else
        set(handles.PlotDataPushbutton,'CData', yellow_image);
        set(handles.PlotPeaksPushbutton,'CData', yellow_image);
        set(handles.DimRedPushbutton,'CData', yellow_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [2 2 2 0 0];
        Processing('error','Incorrect Step',handles)
    end
end
function NeighborsEdit_Callback(hObject, eventdata, handles)
    global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 0 0 0];
    end
    user_entry = round(str2double(get(hObject,'string')));
    if isnan(user_entry)
        set(hObject,'string',1);
    elseif(user_entry>100)
        set(hObject,'string',100);
    elseif(user_entry<=0)
        set(hObject,'string',1);
    else
        set(hObject,'string',num2str(user_entry));
    end
end

%% Clustering Panel
function ClustPopupmenu_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 0 0];
    end
end
function ClusteringPushbutton_Callback(hObject, eventdata, handles)
    global green_image step_execution yellow_image gray_image waiting_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1
        set(handles.ClusteringPushbutton,'CData', waiting_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        Processing('msj','Clustering',handles)
        Experiment=Get_Experiment(handles);
        if isfield(Experiment,'DimensionalReduction')
            Processing('msj','Plotting clustering',handles)
            Experiment=PlotClustering(Experiment,handles);
            Processing('msj','Plotting LLE',handles)
            Experiment=PlotHtoLLE(Experiment,handles);
            Processing('msj','Plotting structure of groups',handles)
            Experiment=PlotGroupStruc(Experiment,handles);
            assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
            set(handles.ClusteringPushbutton,'CData', green_image);            
            step_execution = [1 1 1 1 0];
            Processing('end','Finish Clustering',handles)
        else
            Processing('error','Error Clustering',handles)
        end
        %Processing('end',handles)
    elseif step_execution(1)==1 && step_execution(2)==1 && (step_execution(3)==0 || step_execution(3)==2)
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', yellow_image);
        set(handles.ClusteringPushbutton,'CData', yellow_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 2 2 0];
         Processing('error','Incorrect Step',handles)
    elseif step_execution(1)==1 && (step_execution(2)==0 || step_execution(2)==2) && (step_execution(3)==0 || step_execution(3)==2)
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', yellow_image);
        set(handles.DimRedPushbutton,'CData', yellow_image);
        set(handles.ClusteringPushbutton,'CData', yellow_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 2 2 2 0];
         Processing('error','Incorrect Step',handles)
    else
        set(handles.PlotDataPushbutton,'CData', yellow_image);
        set(handles.PlotPeaksPushbutton,'CData', yellow_image);
        set(handles.DimRedPushbutton,'CData', yellow_image);
        set(handles.ClusteringPushbutton,'CData', yellow_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [2 2 2 2 0];
         Processing('error','Incorrect Step',handles)
    end
end
function LinkagePushbutton_Callback(hObject, eventdata, handles)
    %Processing('start',handles)
    Experiment=Get_Experiment(handles);
    if isfield(Experiment,'DimensionalReduction')
        DataName=Experiment.Data.DataName;
        Sim=Experiment.Peaks.Similarity;
        obj=findobj('name',['Faithfulnes linkage (' DataName ')']);
        if isempty(obj)
            obj=figure('name',['Faithfulnes linkage (' DataName ')'],'numbertitle','off','position',[0 0 600 300]);
        else
            figure(obj);
        end
        clf
        Validation=HBestTree_JP(Sim, obj);
        Experiment.Clustering.Validation=Validation;
        assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
    end
    %Processing('end',handles)
end
function GroupsThSlider_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 0 0];
    end
    slider_value = round(get(handles.GroupsThSlider,'Value'));
    set(handles.GroupsThEdit,'String',slider_value);
end
function GroupsThEdit_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 0 0];
    end   
    user_entry = str2double(get(hObject,'string'));
    if isnan(user_entry)
        Groups=get(handles.GroupsThSlider,'Value');
        set(hObject,'string',num2str(Groups));
        return;
    end
    set(hObject,'string',floor(user_entry));
    user_entry=floor(user_entry);
    if(user_entry>10)
        set(hObject,'string','10');
        set(handles.GroupsThSlider,'Value',10);
    elseif(user_entry<2)
        set(hObject,'string','2');
        set(handles.GroupsThSlider,'Value',2);
    else
        set(handles.GroupsThSlider,'Value',user_entry);
        slider_value = get(handles.GroupsThSlider,'Value');
        set(handles.GroupsThText,'String',strcat('Groups:',32,num2str(slider_value))); 
    end
end
function GroupsPushbutton_Callback(hObject, eventdata, handles)
    %Processing('start',handles)
    Experiment=Get_Experiment(handles);
    if isfield(Experiment,'Clustering')
        Data=Experiment.Peaks.ActiveDataPeaks;
        DataName=Experiment.Data.DataName;
        Adjacency=GetAdjacencyFromPeaks(handles,Data);
        Tree=Experiment.Clustering.Tree;
        obj=findobj('name',['Groups Test (' DataName ')']);
        if isempty(obj)
            obj=figure('name',['Groups Test (' DataName ')'],'numbertitle','off','position',[0 0 600 300]);
        else
            figure(obj);
        end
        clf
        [ClustIdx g]=ClustCorrIdx_JP(Tree,Data,Adjacency,obj);
        
        obj2=findobj('name',['Clustering Coefficient Test (' DataName ')']);
        if isempty(obj2)
            obj2=figure('name',['Clustering Coefficient Test (' DataName ')'],'numbertitle','off','position',[650 70 600 300]);
        else
            figure(obj2);
        end
        clf
        [CoefClustIdx h]=GetClusteringCoef(Tree,'maxclust',Data,Adjacency,obj2);
        
        %
        obj3=findobj('name',['Modularity Test (' DataName ')']);
        if isempty(obj3)
            obj3=figure('name',['Modularity Test (' DataName ')'],'numbertitle','off','position',[650 70 600 300]);
        else
            figure(obj3);
        end
        clf
        [ClustIdx3 g3]=Qidx_JP(Tree,Data,Adjacency,obj3);
        %}
        
        assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
    end
    %Processing('end',handles)
end
function CleanThSlider_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 0 0];
    end
    slider_value = round(get(handles.CleanThSlider,'Value'));
    set(handles.CleanThEdit,'String',slider_value);
end
function CleanThEdit_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 0 0];
    end
    user_entry = str2double(get(hObject,'string'));
    if isnan(user_entry)
        CleanTh=get(handles.CleanThSlider,'Value');
        set(hObject,'string',num2str(CleanTh));
        return;
    end
    user_entry=floor(user_entry);
    if(user_entry>100)
        set(hObject,'string','100');
        set(handles.CleanThSlider,'Value',100);
    elseif(user_entry<50)
        set(hObject,'string','50');
        set(handles.CleanThSlider,'Value',50);
    else
        set(handles.CleanThSlider,'Value',user_entry);
    end
end

%% Network Panel
function CorrMethodPopupmenu_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1 && step_execution(4)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', green_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 1 0];
    end
end
function PlotNetworkPushbutton_Callback(hObject, eventdata, handles)
    global green_image step_execution yellow_image waiting_image 
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1 && step_execution(4)==1
        set(handles.PlotNetworkPushbutton,'CData', waiting_image);
        Processing('msj','Getting experiment',handles)
        Experiment=Get_Experiment(handles);
        if isfield(Experiment,'Structure')
            Processing('msj','Plotting adjacent matrix',handles)
            Experiment=PlotAdjacent(Experiment,handles);
            Processing('msj','Plotting network',handles)
            Experiment=PlotNetwork(Experiment,handles);
            assignin('base',[Experiment.Data.DataName '_Analysis'],Experiment);
            set(handles.PlotNetworkPushbutton,'CData', green_image);
            step_execution = [1 1 1 1 1];
            Processing('end','Finish Plotting Network',handles)
        else
            Processing('error','Plotting Network',handles)
        end
        %Processing('end',handles)
    elseif step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1 && (step_execution(4)==0 || step_execution(4)==2)
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', yellow_image);
        set(handles.PlotNetworkPushbutton,'CData', yellow_image);
        step_execution = [1 1 1 2 2];
         Processing('error','Incorrect Step',handles)
    elseif step_execution(1)==1 && step_execution(2)==1 && (step_execution(3)==0 || step_execution(3)==2) && (step_execution(4)==0 || step_execution(4)==2)
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', yellow_image);
        set(handles.ClusteringPushbutton,'CData', yellow_image);
        set(handles.PlotNetworkPushbutton,'CData', yellow_image);
        step_execution = [1 1 2 2 2];
         Processing('error','Incorrect Step',handles)
    elseif step_execution(1)==1 && (step_execution(2)==0 || step_execution(2)==2) && (step_execution(3)==0 || step_execution(3)==2) && (step_execution(4)==0 || step_execution(4)==2)
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', yellow_image);
        set(handles.DimRedPushbutton,'CData', yellow_image);
        set(handles.ClusteringPushbutton,'CData', yellow_image);
        set(handles.PlotNetworkPushbutton,'CData', yellow_image);
        step_execution = [1 2 2 2 2];
         Processing('error','Incorrect Step',handles)
    else
        set(handles.PlotDataPushbutton,'CData', yellow_image);
        set(handles.PlotPeaksPushbutton,'CData', yellow_image);
        set(handles.DimRedPushbutton,'CData', yellow_image);
        set(handles.ClusteringPushbutton,'CData', yellow_image);
        set(handles.PlotNetworkPushbutton,'CData', yellow_image);
        step_execution = [2 2 2 2 2];
         Processing('error','Incorrect Step',handles)
    end
end
function CorrThSlider_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1 && step_execution(4)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', green_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 1 0];
    end
    slider_value = get(handles.CorrThSlider,'Value');
    if (slider_value)
        slider_value = get(handles.CorrThSlider,'Value');
        set(handles.CorrThEdit,'String',num2str(slider_value*100));
    else
        set(handles.CorrThSlider,'Value',.01);
        set(handles.CorrThEdit,'String',num2str(1));
    end
end
function CorrThEdit_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1 && step_execution(4)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', green_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 1 0];
    end
    user_entry = str2double(get(hObject,'string'));
    if isnan(user_entry)
        CorrTh=get(handles.CorrThSlider,'Value');
        set(hObject,'string',num2str(CorrTh));
        return;
    end
    
    if(user_entry>100)
        set(hObject,'string','100');
        set(handles.CorrThSlider,'Value',1);
    elseif(user_entry<1)
        set(hObject,'string','1');
        set(handles.CorrThSlider,'Value',.01);
    else
        set(handles.CorrThSlider,'Value',user_entry/100);
        set(hObject,'string',num2str(user_entry));
    end
end
function XYSortPopupmenu_Callback(hObject, eventdata, handles)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1 && step_execution(4)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', green_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 1 0];
    end
    %Processing('start',handles)
    %Experiment=Get_Experiment(handles);
    %if isfield(Experiment,'Adjacency')
    %    Experiment=PlotNetwork(Experiment,handles);
    %end
    %Processing('end',handles)
end

%% Colors panel
function GroupColorsPopupmenu_Callback(hObject, eventdata, handles)
    DataPopupmenu_Callback(hObject, eventdata, handles)
end
function XYColorsPopupmenu_Callback(hObject, eventdata, handles)
    DataPopupmenu_Callback(hObject, eventdata, handles)
end
function GetFilePushbutton_Callback(hObject, eventdata, handles)
    global gray_image step_execution
    [FileName,PathName] = uigetfile('*.mat','Select the EXPERIMENT file');
    if isequal(FileName,0)
        %disp('User selected Cancel')
        Processing('end','User selected Cancel',handles)
    else
        set(handles.PlotDataPushbutton,'CData', gray_image);
        set(handles.PlotPeaksPushbutton,'CData', gray_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [0 0 0 0 0];
        
        Processing('msj','Cleaning Workspace',handles)
        %evalin('base', 'clear;')
        Processing('msj','Reading File',handles)
        temp_data=load(fullfile(PathName,FileName)); % will be loaded as structure
        file_variables=fieldnames(temp_data);% get the field names as cell array
        Processing('msj','Loading in Workspace',handles)
        for ii=1:length(file_variables)
            % file_variables{ii} - string of the field name
            % temp_data.(file_variables{ii}) - dynamic field reference
        assignin('base', file_variables{ii}, temp_data.(file_variables{ii}));
        end       
        %Processing('start',handles)
        Processing('msj','Reading Variables',handles)
            DataStrings = evalin('base','who');
            set(handles.DataPopupmenu,'String',[{'---None---'};DataStrings])
            set(handles.DataPopupmenu,'Value',1)
            set(handles.XYPopupmenu,'String',[{'---None---'};DataStrings])
            set(handles.XYPopupmenu,'Value',1)
            set(handles.GroupColorsPopupmenu,'String',[{'R G B Y M C O P Ro Br'};DataStrings])
            set(handles.GroupColorsPopupmenu,'Value',1)
            set(handles.XYColorsPopupmenu,'String',[{'Black'};DataStrings])
            set(handles.XYColorsPopupmenu,'Value',1)
        %Processing('end',handles)
        Processing('end','Loaded File',handles)
    end    
end
% --- Executes when selected object is changed in DimRedUipanel.
function DimRedUipanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in DimRedUipanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', gray_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 0 0 0];
    end
end


% --- Executes when selected object is changed in ClusteringUipanel.
function ClusteringUipanel_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in ClusteringUipanel 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global green_image step_execution gray_image
    if step_execution(1)==1 && step_execution(2)==1 && step_execution(3)==1
        set(handles.PlotDataPushbutton,'CData', green_image);
        set(handles.PlotPeaksPushbutton,'CData', green_image);
        set(handles.DimRedPushbutton,'CData', green_image);
        set(handles.ClusteringPushbutton,'CData', gray_image);
        set(handles.PlotNetworkPushbutton,'CData', gray_image);
        step_execution = [1 1 1 0 0];
    end
end

function WThPushbutton_Callback(hObject, eventdata, handles)

    Experiment=Get_Experiment(handles);
    Data=Experiment.Peaks.ActiveDataPeaks;
    
    N = inputdlg('Enter number of iterations:','Montecarlo test',...
    1,{'1000'});
    N=N{1};
    if isempty(N)
        return;
    else
        WTh=WeightTh(Data, str2num(N), .01);
        set(handles.CorrThEdit,'string',WTh);
        set(handles.CorrThSlider,'value',WTh/100);
        msgbox(['Recommended: ' num2str(WTh)],'Montecarlo test for Correlation Threshold')
    end
end
