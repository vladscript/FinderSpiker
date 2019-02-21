%% Load Tabel **************************************
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end))]; % Finder Spiker Main Folder
% Load File 
[FileName,PathName,MoreFiles] = uigetfile({'*.csv';'*.xlsx'},' Dataset file of ALL Features',...
    'MultiSelect', 'off',CurrentPathOK);
Xraw=readtable([PathName,FileName]);                % Table
Y=categorical(table2array( Xraw(:,1)) );            % Labels
EXPIDs=table2array( Xraw(:,2));                     % Cell of Strings
X=table2array( Xraw(:,3:end) );                     % Dataset
FeatureNames=Xraw.Properties.VariableNames(3:end);  % Feature Names
%% Feature Explorer ********************************

%% Correlation among Features
FeaturesCorr=corr(X);


%% p-values: statistical tests
Nfeat=numel(FeatureNames);
Labels=unique(Y);
Nconditions=numel(Labels);
for f=1:Nfeat
    pMatrix=ones(Nconditions);
    for c=1:Nconditions
        elseCond=setdiff(1:Nconditions,c);
        for e=1:numel(elseCond)
            fprintf('>>Testing %s @ Conditions:\n   %s vs %s p=',FeatureNames{f},char(Labels(c)),char(Labels(elseCond(e))));
            Aindx=find(Y==Labels(c));
            Bindx=find(Y==Labels(elseCond(e)));
            A=X(Aindx,f);
            B=X(Bindx,f);
            [h,p]=ttest2(A,B);
            pMatrix(c,elseCond(e))=p;
            %p=kruskalwallis([A;B],[Y(Aindx);Y(Bindx)])
            fprintf('%3.2f\n',p);
        end
    end
    pause;
end
%% For paired Experiments: Delta Features
% find and compare paired Features

% Calculate Delta-> Define Delta Reference: for Instace: Dyskinesia

