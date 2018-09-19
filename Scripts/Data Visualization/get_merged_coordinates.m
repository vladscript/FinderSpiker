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
function [XY_merged,MetaDataColocaliation]=get_merged_coordinates(Experiment,XY,r)
%% Setup
% Global Variables:
global XY_merged;
% global ColocateIndx;
global MetaDataColocaliation;
Experiment=Experiment(2:end); % delete Slahs from Name ID
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
calciumind= inputdlg('Fluorophoro: ',...
             'Ca++ Indicator', [1 50]);
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
RGBNames={'Red','Green','Blue','Gray (contrast)'};

rgbA=getrgb(ImageAverage,RGBindexesA);
rgbB=getrgb(meanFrame,RGBindexesB);
rgbC=imadd(rgbA,rgbB);

%% Show Images ###################################################
% PLot Figure
getcoord=figure('numbertitle','off',...
    'name',['Colocalization of ',kindcells{1},' cells @ ',Experiment],...
    'Position',[6 151 1099 523]);

%     'keypressfcn',@exit_function,...
%     'ButtonDownFcn',@selec_ROI);
h1=subplot(1,3,1);
imshow(rgbA);
% imcontrast(h1);
h2=subplot(1,3,2);
imshow(rgbB);
h3=subplot(1,3,3);
imshow(rgbC);
% imcontrast(h2);
h1.Position=[-0.34,0.15,1,0.65];
h2.Position=[-0.005,0.15,1,0.65];
h3.Position=[0.33,0.15,1,0.65];
linkaxes([h1,h2,h3],'xy')

title(h1,dyename)
title(h2,calciumind)
title(h3,'Merge & Active Cells')
% Fixing Position
% getcoord.Position=[104 387 1193 283];


% Plot Coordinates & Get ROI pixels
Meshxy_Circle=[]; 
aux1=1;
for n=1:N
%     rectangle('Position',[XY(n,1)-r(n)/2,XY(n,2)-r(n)/2,2*r(n),2*r(n)],...
%             'Curvature',[1 1],'LineWidth',1,'EdgeColor','r');
    disp('Ploting Active Coordinates');
    subplot(h1)
    hold(h1,'on'); % PLOT AT DYE MARKER
    rectangle('Position',[XY(n,1)-r(n),XY(n,2)-r(n),2*r(n),2*r(n)],...
        'Curvature',[1 1],'LineWidth',2,'EdgeColor','k');
    hold(h1,'off'); % PLOT AT DYE MARKER
    
    subplot(h3)
    hold(h3,'on'); % PLOT AT MERGE
    rectangle('Position',[XY(n,1)-r(n),XY(n,2)-r(n),2*r(n),2*r(n)],...
        'Curvature',[1 1],'LineWidth',1,'EdgeColor','r');
    hold(h3,'off'); % PLOT AT MERGE
    
    % ROI pixels
    Mx=XY(n,1)-(r):XY(n,1)+(r); % range in x of square
    My=XY(n,2)-(r):XY(n,2)+(r); % range in y of square
    for i=1:length(Mx)
        % chech if it's in image's limits Xaxis
        if Mx(i)>0 && Mx(i)<=W 
            for j=1: length(My)
                % chech if it's in image's limits Yaxis
                if My(j)>0 && My(j)<=H 
                    % check if it's in circle
                    if (Mx(i)-XY(n,1))^2+(My(j)-XY(n,2))^2<=r(n)^2
                        Meshxy_Circle=[Meshxy_Circle;Mx(i),My(j)];
                        ROIcounter(n)=aux1;
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
'Callback',@ReplotDye,...
'Units','normalized',...
'Value',RGBindexesA,...
'Position',[0.1,0.7,0.1,0.2]);

% Color for Calcium INdicator
rgbChooseCa=uicontrol('Style','popup',...
'String',RGBNames,...
'Callback',@ReplotCa,...
'Value',RGBindexesB,...
'Units','normalized',...
'Position',[0.4,0.7,0.1,0.2]);

% Selector Activator
Selector=uicontrol('Style','togglebutton',...
'String','Select XY',...
'Callback',@getselection,...
'Units','normalized',...
'Position',[0.8 0.85 0.1 0.1]);

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
fhelp=msgbox('Push Select XY and use normal button clicks to add points. A shift-, right-, or double-click adds a final point and ends the selection. Pressing Return or Enter ends the selection without adding a final point. Pressing Backspace or Delete removes the previously selected point.','Help','help');
% IT GETS POINTS ONLY IN THE CIRCLES OF THE COORDINATES
Meshxy_Circle=round(Meshxy_Circle);
X_col=[];Y_col=[];
window_min=0; window_max=255;
% ManageContrast=[];
hb = uicontrol('Style','pushbutton',...
'String','Close & SAVE',...
'Callback',@CloseAndSave);

%% Nested Funtions Magic
% Funtion to Close Figure and Save Stuff 
    function CloseAndSave(~,~)
        % Call global variables
        % global XY_merged;
        % global ColocateIndx;
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
        disp('Rejected Selections: ')
        disp(N_col-numel(ColocateIndx))
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
        SaveImage=frame2im(getframe(hout));
        ImageName=[PathNameVideo,F.Name,'.png'];
        imwrite(SaveImage,ImageName);


        %% Save Colocated Coordinates
        checkname=1;
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
                MetaDataColocaliation.Dye=dyename;
                MetaDataColocaliation.Cells=kindcells;
                save([PathName,FileName],'XY_merged','ColocateIndx','MetaDataColocaliation','-append');
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

% funtion to get RGB image from Grayscale Image:
    function rgbA=getrgb(A,indx)
        rgbA=cat(3,A,A,A);
        if indx==4 % GRAY
            rgbA(:,:,[2,3])=0;
            % rgbA=rgb2gray(rgbA);
        else %  R-G-B
            for k=1:3
                if k~=indx
                    rgbA(:,:,k)=0;
                end
                % Other Combinations: mixes
            end
        end
    end
% Function to Replot FIRST SUBPLOT
    function ReplotDye(~,~)
        RGBindexesA=rgbChooseDye.Value;
        rgbA=getrgb(ImageAverageCopy,RGBindexesA);
        if RGBindexesA<4 
            h1.Children(end).CData=rgbA;
            if and(RGBindexesB<4,RGBindexesA~=RGBindexesB)
                rgbC=imadd(rgbA,rgbB);
                h3.Children(end).CData=rgbC;
            end
        else
            h1.Children(end).CData=ImageAverage;  colormap(gray);
            h1.Children(end).CDataMapping='scaled';
            % uiwait(figsave)
            ManageContrast=imcontrast(h1);
            set(ManageContrast, 'CloseRequestFcn', @getValues);
            waitfor(ManageContrast);
            disp('Dye ok')            
        end
    % Nested in Nested Function to get the valuso of MIN & MAX Gray Histogram
    function getValues(~,~)
        window_min = str2double(get(findobj(ManageContrast, 'tag', 'window min edit'), 'String'));
        window_max = str2double(get(findobj(ManageContrast, 'tag', 'window max edit'), 'String'));
        ImageAverageCopy=ImageAverage;
        ImageAverageCopy(ImageAverage<window_min)=0;
        ImageAverageCopy(ImageAverage>window_max)=255;
        ImageAverageCopy=imadjust(ImageAverageCopy,[window_min,window_max]/255,[]);
        delete(ManageContrast);
        % resume(figsave)
    end  
   end

    function ReplotCa(~,~)
        RGBindexesB=rgbChooseCa.Value;
        rgbB=getrgb(meanFrameCopy,RGBindexesB);
        if RGBindexesB<4
            h2.Children.CData=rgbB;
            if and(RGBindexesA<4,RGBindexesA~=RGBindexesB)
                rgbC=imadd(rgbA,rgbB);
                h3.Children(end).CData=rgbC;
            end
        else
            h2.Children.CData=meanFrame;  colormap(gray);
            h2.Children.CDataMapping='scaled';
            ManageContrast=imcontrast(h2);
            set(ManageContrast, 'CloseRequestFcn', @getValues);
            waitfor(ManageContrast);
            disp('Mean Frame ok')
        end
     % Nested in Nested Function to get the values of MIN & MAX Gray Histogram
     function getValues(~,~)
        window_min = str2double(get(findobj(ManageContrast, 'tag', 'window min edit'), 'String'));
        window_max = str2double(get(findobj(ManageContrast, 'tag', 'window max edit'), 'String'));
        meanFrameCopy=meanFrame;
        meanFrameCopy(meanFrame<window_min)=0;
        meanFrameCopy(meanFrame>window_max)=255;
        meanFrameCopy=imadjust(meanFrameCopy,[window_min,window_max]/255,[]);
        delete(ManageContrast);
        % resume(figsave)
     end          
    end

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
    end

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
            if x_min<1 x_min=1; end
            if y_min<1 y_min=1; end
            if x_max>H x_max=H; end
            if y_max>W y_max=W; end
            window_min = round(0.75*min(min(ImageAverage(x_min:x_max,y_min:y_max))));
            window_max = round(1.25*max(max(ImageAverage(x_min:x_max,y_min:y_max))));
            % DYE IMAGE (neuron marker)
            ImageAverageCopy=ImageAverage;
            ImageAverageCopy=imadjust(ImageAverageCopy,[double(window_min),double(window_max)]/255,[]);
            % REPLOT
            RGBindexesA=rgbChooseDye.Value;
            rgbA=getrgb(ImageAverageCopy,RGBindexesA);
            if RGBindexesA<4 
                h1.Children(end).CData=rgbA;
                if and(RGBindexesB<4,RGBindexesA~=RGBindexesB)
                    rgbC=imadd(rgbA,rgbB);
                    h3.Children(end).CData=rgbC;
                end
            end        
            
            % Calcium Fluorescence (calcium indicator)
            window_min = round(0.75*min(min(meanFrame(x_min:x_max,y_min:y_max))));
            window_max = round(1.25*max(max(meanFrame(x_min:x_max,y_min:y_max))));
            meanFrameCopy=meanFrame;
            meanFrameCopy=imadjust(meanFrameCopy,[double(window_min),double(window_max)]/255,[]);
            % REPLOT
            RGBindexesB=rgbChooseCa.Value;
            rgbB=getrgb(meanFrameCopy,RGBindexesB);
            if RGBindexesA<4 
                h2.Children.CData=rgbB;
                if and(RGBindexesA<4,RGBindexesA~=RGBindexesB)
                    rgbC=imadd(rgbA,rgbB);
                    h3.Children(end).CData=rgbC;
                end
            end        
            fprintf('Cell: %i of %i\n',CellIndx,N);
        end

    end
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
end

%% END OF THE WORLD