% Function to accumulate Rate of Ensmbles(RoE)
% of each Ensmble For each Experiment of the same
% experimental setup, i.e., same experimental conditions

% Input
%   Load experiment manually from different folder: one by one
%   R_Condition:  Selected rasters of the Experiment
%   fs: Sampling Frequency
% Output
%    RoE_ALL:       Rate of Ensembles
%    IEnI:          Inter Ensembles Interval
%    EnD:           Ensembles Duration

%% Setup
% Initial:
runs=1;             % Runs Counter
EXPS={};            % List Of Experiments
% Directory:
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end)),'Processed Data'];
% Load File 
[FileName,PathName,MoreFiles] = uigetfile('*.mat',['Experiment .mat file'],...
    'MultiSelect', 'off',CurrentPathOK);
%% Loop to keep loading files
while MoreFiles
    load([PathName,FileName])
    Nensembles=Features_Condition.Nenscond; % N ensembles in each condition
    NC=numel(Nensembles); % N Conditions
    TheMatrix=Features_Condition.CrossEnsmebleSimm;
    % Start as empty
    if runs==1
        SIM_MATRIX=cell(NC,NC);
    end
    %%  Loop to Accummulate Data
    auxc=1;
    for c=1:NC
        auxd=1;
        Row=Nensembles(c);
        for d=1:NC
            fprintf('Neural Ensembles Simmilarity: %s vs %s\n',Names_Conditions{c},Names_Conditions{d})
            % Get Values
            Col=Nensembles(d);
            MatrixCut=TheMatrix([auxc:auxc+Row-1],[auxd:auxd+Col-1]);
            Vector=MatrixCut(:);
            SIM_MATRIX{c,d}=[SIM_MATRIX{c,d};Vector];
            % pause
            auxd=Col+1;
        end
        auxc=Row+1;
    end
    % Disp Experiments Selected:
    EXPS{runs,1}=Experiment
    CurrentPathOK=PathName;
    runs=runs+1;
    [FileName,PathName,MoreFiles] = uigetfile('*.mat',['Experiment .mat file'],...
    'MultiSelect', 'off',CurrentPathOK);
end
disp('>>end.')
%% Plot & Save Stuff to .mat File at Processed Data Folder
plot_pdf_simmatrix(SIM_MATRIX,Names_Conditions);

okbutton = questdlg('Save data?');
waitfor(okbutton); 
if strcmp('Yes',okbutton)
    % Set Save Name
    timesave=clock;
    TS=num2str(timesave(1:5));
    TS=TS(TS~=' ');
    SaveFile=['\Simm_Matrix_',TS,'.mat'];
    % Select Destiny
    PathSave=uigetdir(CurrentPathOK);
    disp('>>saving data...')
    save([PathSave,SaveFile],'EXPS','Names_Conditions','SIM_MATRIX');
    fprintf('>> Data saved @: %s\n',[PathSave,SaveFile])
else
    fprintf('>>Unsaved data.\n')
end
fprintf('>>Cleaning Workspace: ')
clear
fprintf('done\n')