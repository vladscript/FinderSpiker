%% Get Video Intel: 
% It require to load Experiment 
% W & H must be the same in all experiments
CurrentPath=DefaultPath;
[FileName,PathName] = uigetfile('*.avi',['Select any Video of the Experiment: ',Experiment],...
    'MultiSelect', 'off',CurrentPath);
VideoObj = VideoReader([PathName FileName]);
%F = VideoObj.Duration*VideoObj.FrameRate; % Frames
W = VideoObj.Width;                       % Weight
H = VideoObj.Height;                      % Height

% Select #Video and Condition Video of Experiment

i=1;    % Condition
j=1;    % # Video

X_SPARSE=SIGNALSclean{i,j};
isSIGNAL=isSIGNALS{i,j};

NameOut=[PathName,Names_Conditions{i,j},'_',num2str(j),'_CLEANSED'];

get_clean_video(X_SPARSE,fs,isSIGNAL,XY,r,W,H,NameOut);
