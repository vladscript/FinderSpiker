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
SpacesIndx=find(FileName=='_');
CategoryFeature=FileName(1:SpacesIndx(2)-1);
% Number of Features
[DatasetSize,NFeatures]=size(X);
NameFeatures=Xraw.Properties.VariableNames(3:end);
%% Evaluate Features
% Find Specific Zeros or NaNs
% % RejectData=[];
% % for d=1:DatasetSize
% % %     if ~isempty(X(d,:)==0)
% % %         RejectData=[RejectData;d];
% % %     end
% %     if ~isempty(find(isnan(X(d,:))))    
% %         RejectData=[RejectData;d];
% %     end
% % end
% Set of Features Discardable for each Category
% Raster_Activity
% Ensembles
% Ensembles_de
%% Classifier ALL FEATURES ##########################################
disp('>>Data:')
tabulate(Y)
disp('>>Training...')
Mdl=fitcnb(X,Y,'DistributionNames','kernel');
disp('>>Trained.')
disp('>>Evaluating...')
[Yhat,~]=resubPredict(Mdl);
[C,order]=confusionmat(Y,Yhat)
ECVall=1-numel(find(Y==Yhat))/numel(Y);
fprintf('Cross-Validated Classification Error: %3.1f %%\n',100*ECVall)
disp('>>Evaluated.')
%% Find Failed Experiments:
disp('Missclassified:')
[EXPIDs(Y~=Yhat),Y(Y~=Yhat)]
%% Feature Selection #####################################################
% NFeatSel=2;
Nremove=1;
NFeatSel=NFeatures-Nremove; % Remove-one-Feature Sets
CO=nchoosek(1:NFeatures,NFeatSel); % Sets of Features Indexes
Ncombinations=size(CO,1);
ErrorCV=ones(Ncombinations,1);
for f=1:Ncombinations
    fprintf('Feature Set: %i/%i\n',f,Ncombinations)
    disp('>>Training...')
    Mdl=fitcnb(X(:,CO(f,:)),Y,'DistributionNames','kernel');
    disp('>>Trained.')
    disp('>>Evaluating...')
    [Yhat,~]=resubPredict(Mdl);
    [C,order]=confusionmat(Y,Yhat)
    ECV=1-numel(find(Y==Yhat))/numel(Y);
    fprintf('Cross-Validated Classification Error: %3.1f %%\n',100*ECV);
    ErrorCV(f)=ECV;
    disp('>>Evaluated.')
end
MinError=min(ErrorCV);
BestSetFeatures=find(ErrorCV<=MinError);
RejecetFeat=[];
for n=1:numel(BestSetFeatures)
    RejecetFeat=[RejecetFeat;setdiff(1:NFeatures,CO(BestSetFeatures(n),:))];
end
RejecetFeat=makerowvector(unique(RejecetFeat))';
% Check that the SUBSET Acutally classifies as usin
OKFeatures=setdiff(1:NFeatures,RejecetFeat);
Mdl=fitcnb(X(:,OKFeatures),Y,'DistributionNames','kernel');
disp('>>Trained.')
disp('>>Evaluating...')
[Yhat,~]=resubPredict(Mdl);
ErrorSelFeat=1-numel(find(Y==Yhat))/numel(Y);
fprintf('Cross-Validated Classification Error: %3.1f %%\n',100*ErrorSelFeat);
%% Removing Features
aux=1;
RejectSets={};
TrialError=[];
while ErrorSelFeat<=ECVall
    % ECVall=MinError;
    OKFeatures=setdiff(1:NFeatures,RejecetFeat);
    Nokfeatures=numel(OKFeatures);
    CO=nchoosek(OKFeatures,Nokfeatures-Nremove);
    Ncombinations=size(CO,1);
    ErrorCV=ones(Ncombinations,1);
    for f=1:Ncombinations
        fprintf('Feature Set: %i/%i Trial %i\n',f,Ncombinations,aux)
        disp('>>Training...')
        Mdl=fitcnb(X(:,CO(f,:)),Y,'DistributionNames','kernel');
        disp('>>Trained.')
        disp('>>Evaluating...')
        [Yhat,~]=resubPredict(Mdl);
        [C,order]=confusionmat(Y,Yhat)
        ECV=1-numel(find(Y==Yhat))/numel(Y);
        fprintf('Cross-Validated Classification Error: %3.1f %%\n',100*ECV);
        ErrorCV(f)=ECV;
        disp('>>Evaluated.')
    end
    MinError=min(ErrorCV);
    TrialError(aux)=MinError;
    BestSetFeatures=find(ErrorCV<=MinError);
    % RejecetFeat=[];
    for n=1:numel(BestSetFeatures)
        RejecetFeat=[RejecetFeat;setdiff(1:NFeatures,CO(BestSetFeatures(n),:))'];
    end
    RejecetFeat=unique(RejecetFeat);
    RejectSets{aux}=RejecetFeat;
    OKFeatures=setdiff(1:NFeatures,RejecetFeat);
    Mdl=fitcnb(X(:,OKFeatures),Y,'DistributionNames','kernel');
    disp('>>Trained.')
    disp('>>Evaluating...')
    [Yhat,~]=resubPredict(Mdl);
    [C,order]=confusionmat(Y,Yhat)
    ErrorSelFeat=1-numel(find(Y==Yhat))/numel(Y);
    aux=aux+1;
end
%% Check if it search for less Features
if aux>2
    okaux=aux-2;
elseif aux==1
    RejectSets{1}=RejecetFeat;
    okaux=1;
else
    okaux=1;
end
RejecetFeat=RejectSets{okaux};
OKFeatures=setdiff(1:NFeatures,RejecetFeat);
Mdl=fitcnb(X(:,OKFeatures),Y,'DistributionNames','kernel');
disp('>>Trained.')
disp('>>Evaluating...')
[Yhat,~]=resubPredict(Mdl);
ErrorSelFeat=1-numel(find(Y==Yhat))/numel(Y);
%% Display Exploratory Plot of The Selected Features
fprintf('Cross-Validated Classification Error: %3.1f %%\n',100*ErrorSelFeat);
fprintf('Selected Features: %i of %i:\n',numel(OKFeatures),NFeatures);
NameFeatures(OKFeatures)'
table(order,C)
% pause
%% Boxplots
for n=1:numel(OKFeatures)
    ActualFeature=OKFeatures(n);
    Xdata=X(:,ActualFeature);
    figure;
    boxplot(Xdata,Y)
    title(NameFeatures(OKFeatures(n)))
end
%% Scatter PLots
NFeatComb=nchoosek(OKFeatures,3);
Ncomb=size(NFeatComb,1);
Ncolors=numel(unique(Y));
CM=jet(Ncolors);
Labels=unique(Y);
Colors=zeros(size(Y));
ColorsHatRGB=zeros(numel(Colors),3);
for c=1:Ncolors
    Colors(Y==Labels(c))=c;
    Nrowsid=find(Yhat==Labels(c));
    ColorsHatRGB(Nrowsid,:)=repmat(CM(c,:),numel(Nrowsid),1);
end
% Missclassified * * * *  * * * * 
MissClassData=find(Y~=Yhat);
for n=1:Ncomb
    ComBF=NFeatComb(n,:);
    x=X(:,ComBF(1));
    y=X(:,ComBF(2));
    z=X(:,ComBF(3));
    s=50*ones(size(x));
    figure; 
    hs=scatter3(x,y,z,s,Colors,'filled'); hold on;
    for m=1:numel(MissClassData)
        plot3(x(MissClassData(m)),y(MissClassData(m)),...
            z(MissClassData(m)),...
            'LineStyle','none',...
            'Marker','x',...
            'LineWidth',2,...
            'MarkerSize',15,...
            'MarkerEdgeColor',ColorsHatRGB(MissClassData(m),:)); 
    end
    hs.MarkerEdgeColor='k';
    hold off;
    xlabel(NameFeatures(ComBF(1)));
    ylabel(NameFeatures(ComBF(2)));
    zlabel(NameFeatures(ComBF(3)));
end