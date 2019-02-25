function varargout = Experiment(varargin)
% EXPERIMENT M-file for Experiment.fig
%      EXPERIMENT, by itself, creates a new EXPERIMENT or raises the existing
%      singleton*.
%
%      H = EXPERIMENT returns the handle to a new EXPERIMENT or the handle to
%      the existing singleton*.
%
%      EXPERIMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXPERIMENT.M with the given input arguments.
%
%      EXPERIMENT('Property','Value',...) creates a new EXPERIMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Experiment_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Experiment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help Experiment
% Last Modified by GUIDE v2.5 26-Sep-2012 12:44:07

%% ...::: Initialization :::...
    % Begin initialization code - DO NOT EDIT
    gui_Singleton = 1;
    gui_State = struct('gui_Name',       mfilename, ...
                       'gui_Singleton',  gui_Singleton, ...
                       'gui_OpeningFcn', @Experiment_OpeningFcn, ...
                       'gui_OutputFcn',  @Experiment_OutputFcn, ...
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
% --- Executes just before Experiment is made visible.
function Experiment_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject;
    guidata(hObject, handles);
    if ~isempty(varargin)
        Action=varargin{1};
        switch Action
            case 'Simple'
                %PlotSimple(handles,varargin);
                PlotSimple2(handles,varargin);
            case 'Peaks'
                axes(handles.GroupsAxes)
                plot(rand(5))
        end
    end
end
function varargout = Experiment_OutputFcn(hObject, eventdata, handles) 
    varargout{1} = handles.output;
end

%% Callbacks
%% Toolbar
function SaveUipushtool_ClickedCallback(hObject, eventdata, handles)
    SaveFigure(gcf)
end
%% ...::: Functions :::...

function PlotSimple(handles,options)
% Plot raster
% varargin = {'Simple', Data, DataName, Colors}
    
    % Evaluate varargin
    Data=options{2};
    DataName=options{3};
    Colors=options{4};
    
    [F C]=size(Data);
    [NCellColor NColors]=size(Colors);
    if (NColors==3 && max(max(Colors))<=1 && min(min(Colors))>=0)
        if NCellColor==1
            Colors=repmat(Colors,C,1);
        elseif NCellColor~=C
            Colors=repmat([0 0 0],C,1); % color: black
        end
    else
        Colors=repmat([0 0 0],C,1); % color: black
    end
    % Plot Raster
    axes(handles.RasterAxes); cla
    title(DataName)
    axis([1 F 1 C])
    box
    hold on
    for i=1:C
        idx=find(Data(:,i));
        plot(idx,Data(idx,i)*i,'.','color',Colors(i,:))
    end
    hold off
    set(gca,'XTicklabel','')
    set(gca,'XTick',0)

    % Plot Coactive cells
    axes(handles.CoActiveAxes)
    Co=sum(Data,2);
    plot(Co,'color','black')
    xlim([1 F])

    % Plot Accumulated activity
    axes(handles.ActivityAxes)
    Ac=sum(Data,1);
    plot(Ac,'-','color','black')
    hold on
    for i=1:C
        plot(i,Ac(i),'*','color',Colors(i,:))
    end
    hold off
    set(gca,'XTicklabel','')
    set(gca,'XTick',0)
    xlim([1 C])
    view([90 -90])
end

function PlotSimple2(handles,options)
% Plot raster
% varargin = {'Simple', Data, DataName, Colors}
    
    % Evaluate varargin
    Data=options{2};
    DataName=options{3};
    Colors=options{4};
    
    [F C]=size(Data);
    [NCellColor NColors]=size(Colors);
    if (NColors==3 && max(max(Colors))<=1 && min(min(Colors))>=0)
        if NCellColor==1
            Colors=repmat(Colors,C,1);
        elseif NCellColor~=C
            Colors=repmat([0 0 0],C,1); % color: black
        end
    else
        Colors=repmat([0 0 0],C,1); % color: black
    end
    
    
    % Plot Raster
    obj=findobj('name',['Experiment (' DataName ')']);
    if isempty(obj)
        FigRaster=figure('name',['Experiment (' DataName ')'],'numbertitle','off','position',[0 0 1220 460]);
    else
        FigRaster=figure(obj);
    end
    axes('Tag','RasterAxes','outerposition',[.13 .346 .738 .652])
    title(DataName)
    axis([1 F 1 C]); box; hold on
    for i=1:C
        idx=find(Data(:,i));
        plot(idx,Data(idx,i)*i,'.','color',Colors(i,:))
    end
    hold off
    set(gca,'XTicklabel','','XTick',0)

    % Plot Coactive cells
    axes('Tag','CoactiveAxes')
    set(gca,'outerposition',[.13 0 .738 .326]);
    %set(gca,'title','Coactive neurons');
    Co=sum(Data,2);
    plot(Co,'color','black')
    xlim([1 F])

    % Plot Accumulated activity
    axes('Tag','ActivityAxes')
    set(gca,'outerposition',[.876 .346 .123 .652])
    %set(gca,'title','Activity');
    Ac=sum(Data,1);
    plot(Ac,'-','color','black')
    hold on
    for i=1:C
        plot(i,Ac(i),'*','color',Colors(i,:))
    end
    hold off
    set(gca,'XTicklabel','','XTick',0)
    xlim([1 C])
    view([90 -90])
    
    %axes('Tag','GroupAxes','outerposition',[0 .346 .123 .652]);



end
