%% Script that re-generate Video, Signasl & Raster animation:
% Set Video, Condition and Coordinates
%% Get Video Intel: 
% It require to load Experiment 
% W & H must be the same in all experiments
CurrentPath=pwd;
fprintf('>>Getting Video Frames Size.')
[FileName,PathName] = uigetfile('*.avi',['Select any Video of the Experiment: ',Experiment],...
    'MultiSelect', 'off',CurrentPath);
VideoObj = VideoReader([PathName FileName]);
%F = VideoObj.Duration*VideoObj.FrameRate; % Frames
W = VideoObj.Width;                       % Weight
H = VideoObj.Height;                      % Height

% Select #Video and Condition Video of Experiment

i=1;    % Condition
j=1;    % # Video

% Indexes if the Coordinates ( ALL - XY see @fSIENN)
Indexes=[1,2];

% READ DATA
R=RASTER{j,i};          % Driver
XS=SIGNALSclean{i,j};   % Denoised Sparse Signal
% isSIGNAL=isSIGNALS{i,j};% Signals with Activity
isSIGNAL=find(sum(R,2));% Detected Signals in that Raster
NameOut=[PathName,Names_Conditions{i,j},'_',num2str(j)];
NameOutVidDen=[NameOut,'_VIDCLEANSED'];
%% DENOISED VIDEO 
showmeyourvideo(XS,fs,isSIGNAL,XY,r,W,H,NameOutVidDen);

%% SIGNALS PROCESS


XD=DETSIGNALS{j,i};     % Detrended Signal 
% XS=SIGNALSclean{j,i};   % Denoised Sparse Signal
% isSIGNAL=find(sum(R,2));% Detected Signals in that Raster

% ShowCoordinates to Select;
%   Indexes of XY in -> RASTER & SIGNASL cells
% 
% ShowmeyourCoordinates(XY,r,isSIGNAL);
% c=gcf; 
% XYcellsaxis=gca;
% XYcellsaxis.XLim=[0,W];
% XYcellsaxis.YLim=[0,H];
% fprintf('\nInspect & Select Coordinates of Interest\nAnnotate Indexes in the Script\nThen close the figure\n');
% waitfor(XYcellsaxis);
%% Showtime ###############################################################


IndexesText=[];
for n=1:numel(Indexes)
    IndexesText=[IndexesText,num2str(Indexes(n))];
    if n<numel(Indexes)
        IndexesText=[IndexesText,'_'];
    end
end

% SLIDING WINDOW SIZE: + + + + + + + + + + 
NSECONDS=60;
%   + + + + + + + + + + +  + + + + + + + + 

% SIGNALS
WINDOW_SLIDE=fs*NSECONDS;
% Write Video at Experiment Folder
NameOutSignals=[NameOut,'_Sx_',IndexesText];
showmeyoursignals(XS(Indexes,:),XD(Indexes,:),R(Indexes,:),WINDOW_SLIDE,NameOutSignals)

% RASTER
NameOutRaster=[NameOut,'_Raster'];
showmeyourraster(R,Indexes,fs,NSECONDS,NameOutRaster)
%% END ********************************************************************