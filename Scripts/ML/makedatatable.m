%% Read Data ##############################################################
% X & Y From Multiple Feature Tables
function [X,Y,NamesFeatures,folder_name,MatFileName,EXPIDs]=makedatatable()
X=[]; % Features: Observations x Features
Yall=[]; % Labels
Y=[]; % Labels
readOK = true;
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
Load_Default_Directories;
CurrentPathOK=[Dirpwd(1:slashesindx(end)),FolderNameDataset]; % Finder Spiker Main Folder
ntable=1;
% Table Locations and Label Column Indexes ********************************
while readOK
    % Load File 
    [FileName,CurrentPathOK,MoreFiles] = uigetfile({'*.csv'},' Dataset file of Features',...
        'MultiSelect', 'off',CurrentPathOK);
    if MoreFiles
        Xraw=readtable([CurrentPathOK,FileName]);
        ColumnLabels=Xraw.Properties.VariableNames;
        % Choose label
        [index_var,~] = listdlg('PromptString','Select LABEL Columns name or #',...
                        'SelectionMode','single',...
                        'Name','Label Column',...
                        'ListSize',[300 200],...
                        'ListString',ColumnLabels);
        LabelList(ntable)=index_var;
        Ypart=categorical(table2array( Xraw(:,LabelList(ntable))) );    % Labels
        Yall=[Yall;Ypart];
        % Choose ID column
        [index_var2,index_CHECK2] = listdlg('PromptString','Select ID column CANCEL if no ID',...
                        'SelectionMode','single',...
                        'Name','ID Column',...
                        'ListSize',[600 200],...
                        'ListString',ColumnLabels);
        if index_CHECK2
            IDpart=categorical(table2array( Xraw(:,index_var2)) );
        else
            IDpart=categorical(1:numel(Ypart)); 
            IDpart=IDpart';
        end
        IDsList{ntable}=IDpart;
        NIds(ntable)=numel(IDpart);
        FeatureList{ntable}=setdiff(ColumnLabels,ColumnLabels([index_var,index_var2]));
        TableDirsList{ntable}=[CurrentPathOK,FileName];
        ntable=ntable+1;
        fprintf('>>Table %s  with %i observations\n',FileName,size(Xraw,1));
        LastDirOK=CurrentPathOK;
    else
        readOK=false;
    end
end
NIds=max(NIds);
%% Select Features and Rread/merge tables *********************************

for n=1:ntable-1
    [index_var,~] = listdlg('PromptString','Remove Column Features Cancel if non',...
                        'SelectionMode','multiple',...
                        'Name','Label Column',...
                        'ListSize',[300 200],...
                        'ListString',FeatureList{n});
    % Filter seleted Features
    Features2read{n}=setdiff(FeatureList{n},FeatureList{n}(index_var));
end
%% Homologate labels *******************************************************
LabelsRead=unique(Yall);
name='Input the right labels';
numlines=1;
for n=1:numel(LabelsRead)
    defnames{n}=char(LabelsRead(n));
end
RightLabels=inputdlg(defnames,name,numlines,defnames);

%% Read Table 
fprintf('>>Loading Data< ')
for n=1:ntable-1
    Xraw=readtable(TableDirsList{n});
    Xpart=table2array(Xraw(:,Features2read{n}));           % Features
    Ypart=categorical(table2array( Xraw(:,LabelList(n))) );% Labels
    % Fix Labels
    Ycopy=Ypart;
    for i=1:numel(RightLabels)
        BadLabel=char(defnames(i));
        GoodLabel=char(RightLabels(i));
        Ycopy(Ypart==BadLabel)=GoodLabel;
    end
    % Observations Sorting
    [UniqueIDs{n},Indexes]=unique([IDsList{n},Ycopy],'rows'); % ALREADY SORTED
%     if n>1 && size(Xpart,1)~=NIds
%         % Missing Analyzed Experiments
%         fprintf('>Adding NaNs to missing data\n')
%         Xpart=[Xpart;repmat(NaN,NIds-size(Xpart,1),size(Xpart,2))];
%         strcmp(IDsList{1},IDsList{n})
%     end
    % # DATA #
    X=[X,Xpart(Indexes,:)];
    Y=[Y,Ycopy(Indexes)];
    fprintf(' * ')
end
fprintf('>\n')
%% Check Data
fprintf('\n>>Checking data: ')
for i=1:ntable-1
    for j=i+1:ntable-1
        CHKdataID=setdiff(UniqueIDs{i},UniqueIDs{j},'rows');
        if ~isempty(CHKdataID)
            fprintf('\n****************************************\n')
            fprintf('\n>> CHECK EXPERIMENT IDs !!!!!!\n')
            fprintf('\n****************************************\n')
        else
            fprintf(' + ')
        end
    end
end
fprintf('\n')

%% SUMMARY DATA
% Exp IDs
% Feature Names
[Nobser,Nfeatures]=size(X);
Y = removecats(Y(:,1));
labelconditions=unique(Y);
fprintf('\n +++ Analysis Summary +++ \n')
fprintf('Observations: %i\nFeatures: %i\nLabels: %i ->\n',...
    Nobser,Nfeatures,numel(unique(RightLabels)))
summary(Y(:,1));

%% Save Data
% 
% Save X, Y, Feature Names and EXP_IDs
DataXY=table(UniqueIDs{1}(:,1),Y,X);
NamesFeatures={};
for n=1:ntable-1
    NamesFeatures=[NamesFeatures,Features2read{n}];
end
DataXY.Properties.VariableNames={'EXPIDs','Label','Features'};
EXPIDs=DataXY.EXPIDs;
% Save .mat File

prompt = {'Enter file name to save:'};
dlg_title = 'Save Data';
num_lines = 1;
def = {RightLabels{1}};
FileSaveName= inputdlg(prompt,dlg_title,num_lines,def);
if ~isempty(FileSaveName)
    folder_name = uigetdir(LastDirOK,sprintf( 'Select Folder to save %s.mat Data',FileSaveName{1}));
    folder_name =[folder_name ,'\'];
    MatFileName=[FileSaveName{1},'.mat'];
    save([folder_name,MatFileName],'DataXY','NamesFeatures');
    fprintf('>>Saved X-Y data\n')
else
    fprintf('>>Unsaved X-Y data\n')
end