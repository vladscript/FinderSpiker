% Funtion to Read Images from to find cells that colocates 
% with Calcium Imaging Experiments
% It reads images average and display them from some Dye
% IT also samplea som frames from the CaImaging
% Display a Merge of Both Mean Images: Dye and Calcium
% Save Colocated Coordinates in the .mat File of the Experiment Dataset
% Save Image of Colocated Coordinates of the Merge
% Green Coordinates: Colocated Neurons
% Input
%  Experiment:  Experiment ID:          String ID
%   XY:          Population Coordinates: [x_i,y_i;...]
%   r:           Radious Length:         [r_i,...]
%  Dye Image     Get File using Dialog
%  Calcium Image Get File using Dialog
% 
% Ouput
%               Display Image and Coordinates
% XY_merged     Colocated Coordinates
% MetaDataColocaliation: Dye Marker & Cell Type
% Example:
% Experiment='test'; XY=[100,100]; r=5
% [XY_merged,ColocateIndx]=get_merged_coordinates(Experiment,XY_selected,r);
function [XY_merged,MetaDataColocaliation]=get_merged_coordinates_delux(Experiment,caactind,XY,r)
%% Setup
% Global Variables:
global XY_merged;
% global ColocateIndx;
global MetaDataColocaliation;
Experiment=Experiment(Experiment~='\');     % NAMES PATCH
N=length(XY(:,1));

X_col=[];
Y_col=[];
% ColocateColor=['r','g'];


%% Open Dialogue Box to Read Marker A  & Image: of Colocated Cells :GFP/tdTomato/etc
dyename = inputdlg('Neurons Marker : ',...
             'Colocalizer - DYE', [1 50]);

DefaultPath='C:\Users\Vladimir\Documents\Doctorado\Experimentos';
        if exist(DefaultPath,'dir')==0
            DefaultPath=pwd; % Current Diretory of MATLAB
        end
[FileNameImage,PathNameImage] = uigetfile('*.bmp',[' Pick the ',dyename{1},' Image(s) of ',Experiment],...
            'MultiSelect', 'on',DefaultPath);
% Average Image of Colocated Cells :GFP/tdTomato/etc
if iscell(FileNameImage)
    Ni=numel(FileNameImage);
else
    Ni=1;
end
if Ni>1
    disp('>> Averaging ... ')
    for n=1:Ni
        [A,~]=imread([PathNameImage,FileNameImage{n}]);  % Images
        A=double(A);
        if n==1
            ImageAverage=zeros(size(A));
        end
        ImageAverage=ImageAverage+A;
    end
    disp('>> done.')
    ImageAverage=ImageAverage/Ni;
else
    [ImageAverage,~]=imread([PathNameImage,FileNameImage]);  % Images
end
ImageAverage=uint8(round(ImageAverage));
ImageAverageCopy=ImageAverage;
MetA=whos('ImageAverage');
SizeA=MetA.size;
Nbits=8*MetA.bytes/prod(SizeA);
fprintf(' Image in Grayscale of  %i  bits \n',Nbits);
% [counts,binLocations] = imhist(ImageAverage)
% GrayLimits=find(counts>0);
% MIN_IA=binLocations(GrayLimits(1));
% MAX_IA=binLocations(GrayLimits(end));
% Gray_Map=gray(255);
% B=ind2gray(A,map);                              % Gray Scale Image
[H,W]=size(ImageAverage);                                  % Size of the Image
kindcells = inputdlg('Cell Type : ',...
             'Kind of ... ', [1 50]);
%% Open Dialogue Box to Read Video
% calciumind= inputdlg('Fluorophore: ',...
%              'Ca++ Activity Indicator', [1 50]);

uiwait(msgbox(caactind,'Calcium Imaging Indicator','help'));
calciumind=caactind;

[FileNameVideo,PathNameVideo] = uigetfile('*.avi',[' Pick the ',calciumind{1},' Video(s)',Experiment],...
            'MultiSelect', 'on',PathNameImage);
% Sample Frames of Videos: First, Middle & Final Frame of each Video
if iscell(FileNameVideo)
    Nv=numel(FileNameVideo);
else
    Nv=1;
end
nframes=1;
for n=1:Nv
    if Nv>1
        VideoObj = VideoReader([PathNameVideo,FileNameVideo{n}]);
    else
        VideoObj = VideoReader([PathNameVideo,FileNameVideo]);
    end
    fsvid=VideoObj.FrameRate;       % Frame Rate (phony)
    F = VideoObj.Duration*fsvid;    % Frames
    NS=VideoObj.Duration;           % N Seconds
    Frames2Read=[1/fsvid,round(NS/2),NS-1/fsvid];
    for f=1:length(Frames2Read)
        VideoObj.CurrentTime=Frames2Read(f);
        vidFrame = readFrame(VideoObj);
        vidFrame=double(vidFrame);
        if n==1
            meanFrame=zeros(size(vidFrame));           
        end
        meanFrame=meanFrame+vidFrame;
        disp('Getting (insane in the) Mean Frame... ')
        nframes=nframes+1;
    end
end
meanFrame=meanFrame/(nframes-1);
meanFrame=uint8(round(meanFrame));
meanFrameCopy=meanFrame;
% Output
MetaDataColocaliation.Dye=dyename;
MetaDataColocaliation.Cells=kindcells;

% GSmap=colormap(gray);
% C=ind2gray(meanFrame,GSmap);
% Get Mean Image        
%%% MERGE IMAGES
%Image_Merge=imfuse(ImageAverage,meanFrame,'falsecolor','Scaling','joint', 'ColorChannels', [1,0,2]);   

%% Setup Colors & Names
% Ask for colors RGB -> for marker that colocolize
% RED/GREEN/BLUE->
% Then, ask 1 of the 2
%% MAKE RGB IMAGES & ADD IMAGES-> MERGE
% Load Default Values according to:
% Red     -> 1
% Green   -> 2
% Blue    -> 3
% Default = Green + Red
RGBindexesA=2;
RGBindexesB=1;
RGBNames={'Red','Green','Blue','Gray'};

% Color MAPs
% Nbits=8; % read from image format
Ncolors=2^Nbits;
REDmap=zeros(Ncolors,3);
REDmap(1:Ncolors,1)=[0:Ncolors-1]./Ncolors;

GREENmap=zeros(Ncolors,3);
GREENmap(1:Ncolors,2)=[0:Ncolors-1]./Ncolors;

BLUEmap=zeros(Ncolors,3);
BLUEmap(1:Ncolors,3)=[0:Ncolors-1]./Ncolors;
% Color Matrix Cell
CM{1}=REDmap;
CM{2}=GREENmap;
CM{3}=BLUEmap;
CM{4}=rgb2gray(BLUEmap);
%% Show Images ###################################################
% PLot Figure
getcoord=figure('numbertitle','off',...
    'name',['Colocalization of ',kindcells{1},' cells @ ',Experiment],...
    'Position',[6 151 1099 523],'Color',[0,0,0]);

%     'keypressfcn',@exit_function,...
%     'ButtonDownFcn',@selec_ROI);
h1=subplot(1,3,1);
imshow(ImageAverageCopy);
colormap(h1,CM{RGBindexesA})
% imcontrast(h1);
h2=subplot(1,3,2);
imshow(meanFrameCopy);
colormap(h2,CM{RGBindexesB})

h3=subplot(1,3,3);
rgbA=ind2rgb(ImageAverage,CM{RGBindexesA});
rgbB=ind2rgb(meanFrame,CM{RGBindexesB});
rgbAorigin=ind2rgb(ImageAverage,CM{RGBindexesA});
rgbBorigin=ind2rgb(meanFrame,CM{RGBindexesB});
rgbC=imadd(rgbA,rgbB);
imshow(rgbC);
% imcontrast(h2);

% h1.Position=[-0.34,0.15,W/max([H,W]),H/max([H,W])];
% h2.Position=[-0.005,0.15,W/max([H,W]),H/max([H,W])];
% h3.Position=[0.33,0.15,W/max([H,W]),H/max([H,W])];
h1.Position=[0.0,0.15,0.3,0.5];
h2.Position=[0.3,0.15,0.3,0.5];
h3.Position=[0.6,0.15,0.3,0.5];

linkaxes([h1,h2,h3],'xy')

title(h1,dyename,'Color',[0.7,0.7,0.9])
title(h2,calciumind,'Color',[0.7,0.7,0.9])
title(h3,'Merge & Active Cells','Color',[0.7,0.7,0.9])
% Fixing Position
% getcoord.Position=[104 387 1193 283];


% Plot Coordinates & Get ROI pixels
Meshxy_Circle=[]; 
aux1=1;
for nc=1:N
%     rectangle('Position',[XY(n,1)-r(n)/2,XY(n,2)-r(n)/2,2*r(n),2*r(n)],...
%             'Curvature',[1 1],'LineWidth',1,'EdgeColor','r');
    disp('Ploting Active Coordinates');
    subplot(h1)
    hold(h1,'on'); % PLOT AT DYE MARKER
    rectangle('Position',[XY(nc,1)-r(nc),XY(nc,2)-r(nc),2*r(nc),2*r(nc)],...
        'Curvature',[1 1],'LineWidth',2,'EdgeColor','k');
    hold(h1,'off'); % PLOT AT DYE MARKER
    
    subplot(h3)
    hold(h3,'on'); % PLOT AT MERGE
    rectangle('Position',[XY(nc,1)-r(nc),XY(nc,2)-r(nc),2*r(nc),2*r(nc)],...
        'Curvature',[1 1],'LineWidth',1,'EdgeColor','r');
    hold(h3,'off'); % PLOT AT MERGE
    
    % ROI pixels
    Mx=XY(nc,1)-(r):XY(nc,1)+(r); % range in x of square
    My=XY(nc,2)-(r):XY(nc,2)+(r); % range in y of square
    for i=1:length(Mx)
        % chech if it's in image's limits Xaxis
        if Mx(i)>0 && Mx(i)<=W 
            for j=1: length(My)
                % chech if it's in image's limits Yaxis
                if My(j)>0 && My(j)<=H 
                    % check if it's in circle
                    if (Mx(i)-XY(nc,1))^2+(My(j)-XY(nc,2))^2<=r(nc)^2
                        Meshxy_Circle=[Meshxy_Circle;Mx(i),My(j)];
                        ROIcounter(nc)=aux1;
                        aux1=aux1+1;
                    end
                end
            end
        end
    end
end
hold(h3,'on'); % PLOT AT MERGE: stays ON to put MERGED ONES!!!!!
%% SET GUI Handles
% Color for Dye
rgbChooseDye=uicontrol('Style','popup',...
'String',RGBNames,...
'Callback',{@ReplotImage,'Dye'},...
'Units','normalized',...
'Value',RGBindexesA,...
'Position',[0.1,0.7,0.1,0.2]);

% Imcontrast for Dye Image
DyeContrast=uicontrol('Style','togglebutton',...
'String','Enhance',...
'Units','normalized',...
'Callback',{@InitializeImContrast,'Dye'},...
'Position',[0.3 0.9 0.075 0.1]);


% Color for Calcium INdicator
rgbChooseCa=uicontrol('Style','popup',...
'String',RGBNames,...
'Callback',{@ReplotImage,'Ca'},...
'Value',RGBindexesB,...
'Units','normalized',...
'Position',[0.4,0.7,0.1,0.2]);

% Imcontrast for Dye Image
FluoContrast=uicontrol('Style','togglebutton',...
'String','Enhance',...
'Callback',{@InitializeImContrast,'Ca'},...
'Units','normalized',...
'Position',[0.6 0.9 0.095 0.1]);

% Selector Activator
Selector=uicontrol('Style','togglebutton',...
'String','Select XY',...
'Callback',@getselection,...
'Units','normalized',...
'Position',[0.9 0.9 0.1 0.1]);

% Cell Navigator
CellNavigator=uicontrol('Style','slider',...
'Min',0,'Max',N,...
'Value',0,...
'SliderStep',[1/(N),5/(N)],...
'Callback',@cellnavigation,...
'Units','normalized',...
'Position',[0.35,0.05,0.3,0.05]);


% WAIT TO SEE IMAGE
% pause;
%% Get Points 
% Help Message
fhelp=msgbox({'Push Select XY and left clic to add points in the Merge Panel',...
    '',...
    'Pressing Supr or Delete to remove selected point.',...
    '',...
    'Press Enter to add coordinate',...
    '',...
    'Navigate each Coordinate with the Scrollbar below',...
    '',...
    'Press Enhace to adjust contrast.',...
    '',...
    'Press Close & Save to update changes on the .mat file'},'How to','help');
% IT GETS POINTS ONLY IN THE CIRCLES OF THE COORDINATES
Meshxy_Circle=round(Meshxy_Circle);
X_col=[];Y_col=[];
window_min=0; window_max=255;
% ManageContrast=[];
hb = uicontrol('Style','pushbutton',...
'String','Close & SAVE',...
'Units','pixels',...
'Position',[20,20,100,20],...
'Callback',@CloseAndSave);

%% Nested Funtions Magic
% Funtion to Close Figure and Save Stuff 


% funtion to get RGB image from Grayscale Image:
%     function rgbA=getrgb(A,indx)
%         rgbA=cat(3,A,A,A);
%         if indx==4 % GRAY
%             rgbA(:,:,[2,3])=0;
%             % rgbA=rgb2gray(rgbA);
%         else %  R-G-B
%             for k=1:3
%                 if k~=indx
%                     rgbA(:,:,k)=0;
%                 end
%                 % Other Combinations: mixes
%             end
%         end
%     end
    %% Function to Replot FIRST SUBPLOT


%     function ReplotCa(~,~)
%         RGBindexesB=rgbChooseCa.Value;
%         rgbB=getrgb(meanFrameCopy,RGBindexesB);
%         if RGBindexesB<4
%             h2.Children.CData=rgbB;
%             if and(RGBindexesA<4,RGBindexesA~=RGBindexesB)
%                 rgbC=imadd(rgbA,rgbB);
%                 h3.Children(end).CData=rgbC;
%             end
%         else
%             h2.Children.CData=meanFrame;  colormap(gray);
%             h2.Children.CDataMapping='scaled';
%             ManageContrast=imcontrast(h2);
%             set(ManageContrast, 'CloseRequestFcn', @getValues);
%             waitfor(ManageContrast);
%             disp('Mean Frame ok')
%         end
%      % Nested in Nested Function to get the values of MIN & MAX Gray Histogram
%      function getValues(~,~)
%         window_min = str2double(get(findobj(ManageContrast, 'tag', 'window min edit'), 'String'));
%         window_max = str2double(get(findobj(ManageContrast, 'tag', 'window max edit'), 'String'));
%         meanFrameCopy=meanFrame;
%         meanFrameCopy(meanFrame<window_min)=0;
%         meanFrameCopy(meanFrame>window_max)=255;
%         meanFrameCopy=imadjust(meanFrameCopy,[window_min,window_max]/255,[]);
%         delete(ManageContrast);
%         % resume(figsave)
%      end          
%     end




%     function exit_function(getcoord,~,~)
%        key=get(getcoord,'CurrentKey');
%        if strcmp(key,'e') % Exit
%             close(getcoord);
%        end
%     end
 
%     function selec_ROI(~,~)
%         [X_col,Y_col]=getpts;
%         analyze_colocalization();
%     end

%     function analyze_colocalization()
%         N_col=length(X_col);
%         for i=1:N_col
%             for j=1:N
%                 dis(j)=abs(X_col(i)-XY(j,1))+abs(Y_col(i)-XY(j,2));
%             end
%             [~,locind]=min(dis);
%             ColocateIndx(locind)=1;
%             XY_merged=[XY_merged;XY(locind,:)];
%         end
%     end

%% NESTED FUNTIONS
%  REPLOT DYE IMAGE ****************************************
function ReplotImage(~,~,axisIndx)
    % Setup
    if strcmp(axisIndx,'Dye')
        disp('>>Dye')
        axish=h1;
        ImageOriginal=ImageAverage;
        % Output    
        RGBindexes=rgbChooseDye.Value;
        RGBindexesAlter=rgbChooseCa.Value;
        RGBindexesA=RGBindexes;
        RGBindexesB=RGBindexesAlter;
    else
        disp('>>Calcium Fluorescence');
        ImageOriginal=meanFrame;
        RGBindexes=rgbChooseCa.Value;
        RGBindexesAlter=rgbChooseDye.Value;
        axish=h2;
        % Output
        RGBindexesB=RGBindexes;
        RGBindexesA=RGBindexesAlter;
    end
    colormap(axish,CM{RGBindexes});
    % RE-SET ORIGNAL rgb IMAGES
    rgbAorigin=ind2rgb(ImageAverage,CM{RGBindexesA});
    rgbBorigin=ind2rgb(meanFrame,CM{RGBindexesB});
    if RGBindexes~=RGBindexesAlter
        rgbA=ind2rgb(ImageAverageCopy,CM{RGBindexesA});
        rgbB=ind2rgb(meanFrameCopy,CM{RGBindexesB});
        rgbC=imadd(rgbA,rgbB);
        h3.Children(end).CData=rgbC;
        h3.Children(end).CDataMapping='scaled';
    end

    
    % rgbA=getrgb(ImageAverageCopy,RGBindexesA);
%     if RGBindexes<numel(RGBNames)
%         colormap(axish,CM{RGBindexes})
%         % h1.Children(end).CData=rgbA;
%         if and(RGBindexesAlter<numel(RGBNames),RGBindexes~=RGBindexesAlter)
%             rgbA=ind2rgb(ImageAverageCopy,CM{RGBindexes});
%             rgbB=ind2rgb(meanFrameCopy,CM{RGBindexesAlter});
%             rgbC=imadd(rgbA,rgbB);
%             h3.Children(end).CData=rgbC;
%             h3.Children(end).CDataMapping='scaled';
%         end
%     else
%         axish.Children(end).CData=ImageOriginal;  colormap(gray);
%         axish.Children(end).CDataMapping='scaled';
%         % uiwait(figsave)
%         % ManageContrast=imcontrast(h1);
%         % set(ManageContrast, 'CloseRequestFcn', @getValues);
%         % waitfor(ManageContrast);
%         disp('Gray Matter(s)')            
%     end
end

% CONTRAST *************************************************
function InitializeImContrast(~,~,OriginContrast)
    disp('>> Adjust Contrast of ')
    if strcmp(OriginContrast,'Dye')
        disp('>> Dye')
        axish=h1;
        % Restart Original Contrast
        ImageAverageCopy=ImageAverage;
        ImageOriginal=ImageAverage;
        ImageCopy=ImageAverageCopy;
    else
        disp('>> Calcium Activity Fluorescence')
        axish=h2;
        % Restart Original Contrast
        meanFrameCopy=meanFrame;
        ImageOriginal=meanFrame;
        ImageCopy=meanFrameCopy;
    end
    % Re Plot Original Stuff
    axish.Children(end).CData=ImageOriginal;
    axish.Children(end).CDataMapping='scaled';
    % Restart Contrast
    ReplotImage(axish,1,OriginContrast)
    % Adjust CONTRAST GUI
    disp('>>Initializing Image Contrast Adjust...')
    ManageContrast=imcontrast(axish);
    ManageContrast.Position=[ 948,287,423,350];
    set(ManageContrast, 'CloseRequestFcn', @getValues);
    waitfor(ManageContrast);
    disp('...>>Contrast Adjusted.')
    % Replot Stuff and Region
    ReplotImage(axish,1,OriginContrast)
    % Nested in Nested Function to get the values of MIN & MAX Gray Histogram
    function getValues(~,~)
        disp('>> Updating Min-Max Contrast Values')
        window_min = str2double(get(findobj(ManageContrast, 'tag', 'window min edit'), 'String'));
        window_max = str2double(get(findobj(ManageContrast, 'tag', 'window max edit'), 'String'));
        ImageCopy=ImageOriginal;
        ImageCopy(ImageOriginal<window_min)=0;
        ImageCopy(ImageOriginal>window_max)=255;
        ImageCopy=imadjust(ImageCopy,[window_min,window_max]/255,[]);
        delete(ManageContrast);
        % Output
        if strcmp(OriginContrast,'Dye')
            disp('>> Dye Updated')
            ImageAverageCopy=ImageCopy;
        else
            disp('>> Frame Updated')
            meanFrameCopy=ImageCopy;
        end
        disp('>> Done')
        % resume(figsave)
    end    
end

% CELL NAVIGATION ******************************************
function cellnavigation(CellNavigator,~)
    CellIndx=round(CellNavigator.Value);
    if CellIndx>0 && CellIndx<=N
        x_min=XY(CellIndx,1)-3*r(CellIndx)+1;
        x_max=XY(CellIndx,1)+3*r(CellIndx);
        y_min=XY(CellIndx,2)-3*r(CellIndx)+1;
        y_max=XY(CellIndx,2)+3*r(CellIndx);
        axis([x_min,x_max,y_min,y_max])
        % For DYE h1 ##################################################
        % ONLY ROI - - - - - - - - - - - - - - - - - - - - - - - - - - 
        x_min=x_min+r(CellIndx);
        y_min=y_min+r(CellIndx);
        x_max=x_max-r(CellIndx);
        y_max=y_max-r(CellIndx);
        % Take care of negative or BIGGER pixels:!!!!!!!
        if x_min<1; x_min=1; end
        if y_min<1; y_min=1; end
        if x_max>W; x_max=W; end
        if y_max>H; y_max=H; end
        % Adjust Values
        disp('>>Adjusting Contrast ... ')
        
        % DYE IMAGE (neuron marker) #######################################
        AuxDye=rgb2ind(rgbAorigin,CM{RGBindexesA});
        [window_min,window_max]=get_limits(AuxDye(y_min:y_max,x_min:x_max));
        ImageAverageCopy=imadjust(AuxDye,[window_min,window_max]/255,[]);
        % CA++ Indicator IMAGE ############################################
        AuxFluo=rgb2ind(rgbBorigin,CM{RGBindexesB});
        [window_min,window_max]=get_limits(AuxFluo(y_min:y_max,x_min:x_max));
        meanFrameCopy=imadjust(AuxFluo,[window_min,window_max]/255,[]);
        
%         ImageAverageCopy=imadjust(ImageAverageCopy,[window_min,window_max]/255,[]);
        % REPLOT
        % RGBindexesA=rgbChooseDye.Value;
        % rgbA=getrgb(ImageAverageCopy,RGBindexesA);
        rgbA=ind2rgb(ImageAverageCopy,CM{RGBindexesA});
        rgbB=ind2rgb(meanFrameCopy,CM{RGBindexesB});
        if RGBindexesA<numel(RGBNames)
            h1.Children(end).CData=ImageAverageCopy;
            if and(RGBindexesB<4,RGBindexesA~=RGBindexesB)
                rgbC=imadd(rgbA,rgbB);
                h3.Children(end).CData=rgbC;
            end
        end        

        % Calcium Fluorescence (calcium indicator)
        
        % [window_min,window_max]=get_limits(meanFrame(y_min:y_max,x_min:x_max));
        
        % REPLOT
        RGBindexesB=rgbChooseCa.Value;
        rgbA=ind2rgb(ImageAverageCopy,CM{RGBindexesA});
        rgbB=ind2rgb(meanFrameCopy,CM{RGBindexesB});
        % rgbB=getrgb(meanFrameCopy,RGBindexesB);
        if RGBindexesA<4 
            h2.Children.CData=meanFrameCopy;
            if and(RGBindexesA<4,RGBindexesA~=RGBindexesB)
                rgbC=imadd(rgbA,rgbB);
                h3.Children(end).CData=rgbC;
            end
        end        
        fprintf('Cell: %i of %i\n',CellIndx,N);
    end
    % Nested Inception Function
    function[minlimit,maxlimit]=get_limits(ImageSegment)
        disp('>>Recalculating Limits')
%         [Iind,mapind]=gray2ind(ImageSegment,2^Nbits);
        grayvalues=double(ImageSegment(:));
        if numel(grayvalues)>4
            [px,binx]=ksdensity(grayvalues,linspace(min(grayvalues),max(grayvalues),100));
            [~,MaxGray,WidthGray]=findpeaks(px,binx,'Npeaks',1,'SortStr','descend');
            if ~isempty(MaxGray)
                minlimit=MaxGray-1.25*WidthGray;
                maxlimit=MaxGray+1.25*WidthGray;
            else
                minlimit=min(grayvalues);
                maxlimit=max(grayvalues);
            end
        else
            minlimit=min(grayvalues);
            maxlimit=max(grayvalues);
        end
        disp('>>Done')
    end
end

% POINTING COORDINATES
function getselection(~,~)
    if Selector.Value==1
       % Accumulate Coordinates:
        [x_col,y_col]=getpts(getcoord);
        X_col=[X_col;x_col];
        Y_col=[Y_col;y_col];
        Selector.Value=0;
        % set a marker for selected coordinates
        for ni=1:numel(x_col)
            plot(h3,x_col(ni),y_col(ni),'*w')
        end
    end
    % Maybe Save Here Captures of Images
end

% CLOSE AND SAVE *******************************************
function CloseAndSave(~,~)
    % Call global variables
    % global XY_merged;
    % global ColocateIndx;
    axis tight; msgremind=msgbox('Adjust Contrast to Save Image');
    waitfor(msgremind);
    disp('Selection Ready...')
    % [X_col,Y_col]=getpts(getcoord);
    X_col=round(X_col);
    Y_col=round(Y_col);
    % CHeck POint are in ROIs -----------------------------------------------
    N_col=length(X_col);
    AcceptedIndx=[];
    for k=1:N_col
        disp('Checking Point:')
        [~,location]=ismember([X_col(k),Y_col(k)],Meshxy_Circle,'rows');
        % if ismember([X_col(i),Y_col(i)],Meshxy_Circle,'rows')
        % if ~isempty(location)
        if ~location==0
            AuxInd=find(ROIcounter>location);
            % ColocateIndx(AuxInd(1))=1;
            XY_merged=[XY_merged;XY(AuxInd(1),:)];
            AcceptedIndx=[AcceptedIndx,k];
            disp('... Accepted')
        else
            disp('... Rejected')
        end
    end
    % Ignore Repeated Pints in the same ROI
    [XY_merged,~,~]=unique(XY_merged,'rows','stable');
    % ColocateIndx=ColocateIndx(IndOK);
    % disp('Rejected Selections: ')
    % disp(N_col-numel(ColocateIndx))
    close(getcoord);
    % Final Output:
    F=figure;
    if isempty(AcceptedIndx)
        disp('***********NO CO-LOCALIZATION***************')
        hold off;
        imshow(rgbC);
    else
        hold off;
        imshow(rgbC);
        % imcontrast;
        hold on;
        % Plot Coordinates & Get ROI pixels
        for n=1:length(AcceptedIndx)
            rectangle('Position',[XY_merged(n,1)-r(n),XY_merged(n,2)-r(n),2*r(n),2*r(n)],...
            'Curvature',[1 1],'LineWidth',1,'EdgeColor','g');
        end
        % Save IMAGE Manually
    end
    F.Name=[kindcells{1},'_in_',Experiment];
    hout=gca;
    hout.Position=[0,0,1,1];
    disp('Saving Merge Image at:')
    SaveImage=frame2im(getframe(hout));
    ImageName=[PathNameVideo,F.Name,'.png'];
    disp(ImageName)
    imwrite(SaveImage,ImageName);
    disp('Saved.')

    %% Save Colocated Coordinates (INACTIVE)
    checkname=0;
    while checkname==1
        DP=pwd;
        Slashes=find(DP=='\');
        DefaultPath=[DP(1:Slashes(end)),'Processed Data'];
        if exist(DefaultPath,'dir')==0
            DefaultPath=pwd; % Current Diretory of MATLAB
        end
        [FileName,PathName] = uigetfile('*.mat',[' Pick the Analysis File ',Experiment],...
            'MultiSelect', 'off',DefaultPath);
        dotindex=find(FileName=='.');
        if strcmp(FileName(1:dotindex-1),Experiment(1:end))
            checkname=0;
            % SAVE DATA
            
            save([PathName,FileName],'XY_merged','MetaDataColocaliation','-append');
            disp([Experiment,'   -> UPDATED (Merged Coordinates)'])
        elseif FileName==0
            checkname=0;
            disp('....CANCELLED')
        else
            disp('Not the same Experiment!')
            disp('Try again!')
        end
    end    
end
end %% END OF THE WORLD