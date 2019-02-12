% Script To Evaluate Which Category of Features
% Classify the Best the Measured Conditions
%% Setup
%% Load Dataset
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end))]; % Finder Spiker Main Folder
% Load File 
[FileName,PathName,MoreFiles] = uigetfile('*.csv',' Dataset file',...
    'MultiSelect', 'off',CurrentPathOK);

Xraw=readtable([PathName,FileName]);
% Labels
Y=categorical(table2array( Xraw(:,1)) );    % Categorical ARRAYS
EXPIDs=table2array( Xraw(:,2));             % Cell of Strings
X=table2array( Xraw(:,3:end) );
%% Get Intel of the Feature Category
% Category
% Number of Features
%% Evaluate Features
% Find Specific Zeros or NaNs
%% Classifier
disp('>>Data:')
tabulate(Y)
disp('>>Training...')
Mdl=fitcnb(X,Y,'DistributionNames','kernel');
disp('>>Trained.')
disp('>>Evaluating...')
[Yhat,~]=resubPredict(Mdl);
[C,order]=confusionmat(Y,Yhat)
ECV=1-numel(find(Y==Yhat))/numel(Y);
fprintf('Cross-Validated Classification Error: %3.1f %%\n',100*ECV)
disp('>>Evaluated.')
%% Find Failed Experiments:
disp('Missclassified:')
[EXPIDs(Y~=Yhat),Y(Y~=Yhat)]