%% PCA & SVM Classification ###############################################
% 
% Read Table(s) and Select Features from Datasets:
%  - Raster_Activity_Dataset_XXXX.csv
%  - General_Ensembles_Dataset_YYYY.CSV
%  - Network_Features_Dataset_ZZZZ.csv
%  - Or CSV tables with Features/ Columns and ID columns
% IMPORTANT:
% It assumes IDs are the same among tables (for sorting purposes)
% 
%% Setup:
clc; clear; close all;
Import_FinderSpiker;
Nsim=10; % Repetitions for SVM

%% Select Load Data mode #################################################
% Select kind of file:
choiceFile=0;
fprintf('>>Select data file or make data table\n');
fprintf('\n 1)Load Data File:\n')
fprintf('   Select .mat file with already selected features\n')
fprintf('\n 2)Make Data Table_\n')
fprintf('   Select Features from CSV files\n')
while choiceFile==0
    choiceFile = menu('Choose a Data File or Build Data Table','Load Data File','Make Data File');
end
switch choiceFile
    case 1
        fprintf('>> Load Data File \n')
        % Load .mat File
        Dirpwd=pwd;
        slashesindx=find(Dirpwd=='\');
        CurrentPathOK=[Dirpwd(1:slashesindx(end))];
        [FileMat,FolderName]=uigetfile(CurrentPathOK);
        load([FolderName,FileMat])
        X=DataXY.Features;
        Y=DataXY.Label;
        clear DataXY;
    case 2
        fprintf('>>Make Data Table\n')
        % Run Interface to load tables and select Features
        [X,Y,NamesFeatures,FolderName,FileMat]=makedatatable();
end

%% Color for Labels #######################################################
[Nobser,Nfeatures]=size(X);
labelconditions=unique(Y);

for n=1:numel(labelconditions)
    condition_labels{n}=char(labelconditions(n));
end

[CM,ColorIndx]=Color_Selector(condition_labels);
% 
%% PCA ####################################################################
% Each principal component is a linear combination of the original variables.
% All the principal components are orthogonal to each other, so there is no redundant information.
% But it is commonplace for the sum of the variances of the first few principal components to exceed 80% of the total variance of the original data.
% Input data for which to compute the principal components, 
% specified as an n-by-p matrix. 
% SVD: Single Value Decomposition
% Rows of X correspond to observations and columns to variables.
% The variable weights are the inverse of sample variance.
% 
% score, contains the coordinates of the original data in the new coordinate system defined by the principal components
% wcoefs, contains the coefficients of the principal components

fprintf('>>PCA : ... ')
[wcoefs,score,latent,~,explained] = pca(X,'algorithm','svd',...
                        'Centered','off',...
                        'VariableWeights','variance');
fprintf('done.\n')
% Orthoganility:
fprintf('>>Checking Orthonormality: ... ')
coefforth = inv(diag(std(X)))*wcoefs;
c3 = coefforth(:,1:3);
I = c3'*c3
fprintf('done.\n')
% % Display 2 PCA components scatter plot:
% figure
% gscatter(score(:,1),score(:,2),Y,CM(ColorIndx,:),'dso',10);
% xlabel('1st Principal Component')
% ylabel('2nd Principal Component')
% % Rows of score correspond to observations, 
% and columns to components.
ThrehsolVariance=0.95; % *100=% of variance explained by PCA
VarExplained=cumsum(latent)./sum(latent);
Npcs=find(VarExplained>=ThrehsolVariance,1);
if Npcs==1 Npcs=2; end
% Only the first 95% of the cumulative distribution is displayed
figure; 
pareto(explained); hold on
xlabel('Principal Component')
ylabel('Variance Explained (%)')
% plot([Npcs,Npcs],[0,100],'r');

Xpca=score(:,1:Npcs);

%% Multiclass Classification #############################################
% Features:     Xpca,   dims:   Nobservations x Npcs
% Labels:       Y,      dims:   Nobservations x 1
fprintf('\n>> For every kernel and pair of PCs, calculate: ')
fprintf('\n>> LOOCV (leave-one-out cross validation) Accuracy:')
% Test 3 function kernels: linear, parabolic and gaussian
kernels2test={'linear';'gaussian';'polynomial'};
Combo=combnk(1:Npcs,2);
Ncombo=size(Combo,1);
ALLPCAs=[Combo,zeros(Ncombo,3)];
% SVM parameters:
% If you specify KernelScale 'auto', then the software uses a heuristic 
% procedure to select the scale value. The heuristic procedure uses subsampling. 
% Therefore, to reproduce results, set a random number seed using rng before 
% training the classifier.
for k=1:numel(kernels2test)
    % generate SVM template
    fprintf('\n>>Making SVM Template with kernel: %s ',kernels2test{k})
    if k<3
        template = templateSVM('KernelFunction',kernels2test{k}, 'PolynomialOrder', [],...
        'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);    
%         'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);    
        
    else
        template = templateSVM('KernelFunction',kernels2test{k}, 'PolynomialOrder', 2,...
        'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);    
%         'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);
        fprintf(' of order : 2')
    end
    fprintf('\n')
    % SVM Classificator using 2 PCAs: all combinations
    for c=1:Ncombo
        i=ALLPCAs(c,1); j=ALLPCAs(c,2);
        % Make Cross Validation
        Xdata=Xpca(:,[i,j]);
        Yhatsvm=categorical;
        fprintf('\n>>Model %i/%i,kernel %s, PCs: %i vs %i -> \n',c,Ncombo,kernels2test{k},i,j);
        for ns=1:Nsim
            Mdl = fitcecoc(Xdata,Y, 'Learners', template, 'Coding',...
                    'binarycomplete', 'PredictorNames', {['PCA_',num2str(i)],['PCA_',num2str(j)]}, ...
                    'ResponseName', 'Y','FitPosterior',1, ...
                    'ClassNames', labelconditions);
            Yhatsvm = predict(Mdl,Xdata);
            fprintf('*')
            validationAccuracyOK(ns)=sum(Yhatsvm==Y)/numel(Y);
        end
        fprintf('\n')
        LOOCVs{c,k}=validationAccuracyOK;
        ALLPCAs(c,2+k)=mean(validationAccuracyOK);
        fprintf('\n>>Mean Accuracy: %3.2f\n',mean(validationAccuracyOK));        
    end
end
% fprintf(repmat('*',10,1)); fprintf('[complete]\n');
% GET BEST Kernel multiclassificator
if size(ALLPCAs(:,3:end),1)>1
    MaXeachKernel=max(ALLPCAs(:,3:end));
else
    MaXeachKernel=ALLPCAs(:,3:end);
end
[AccuracySVM,kernelindx]=max(MaXeachKernel);

if kernelindx<3
    template = templateSVM('KernelFunction',kernels2test{kernelindx}, 'PolynomialOrder', [],...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);    
%     'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);    
else
    template = templateSVM('KernelFunction',kernels2test{kernelindx}, 'PolynomialOrder', 2,...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);    
%     'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);
end
% BEST PCAs *************************
[~,nModel]=max(ALLPCAs(:,2+kernelindx));
columns2boxplot(LOOCVs{nModel,1}',LOOCVs{nModel,2}',LOOCVs{nModel,3}',{'Lineal','Gaussian','Polynomial'})
ONEpca=ALLPCAs(nModel(1),1);
TWOpca=ALLPCAs(nModel(1),2);
% LOOCV *************************
% clear Yhatsvm
fprintf('Model| Leave-One-Out Cross Validations \n')
for ns=1:Nsim
    Ypre=categorical;
    fprintf('%i/%i ',ns,Nsim)
    for n=1:Nobser
        Xtrain=Xpca(setdiff(1:Nobser,n),[ONEpca,TWOpca]);
        Ytrain=Y(setdiff(1:Nobser,n),:);
        Xtest=Xpca(n,[ONEpca,TWOpca]);
        Ytest=Y(n,:);
        Yhatsvm=categorical;
        Mdl = fitcecoc(Xtrain,Ytrain, 'Learners', template, 'Coding',...
            'binarycomplete', 'PredictorNames', {['PCA_',num2str(i)],['PCA_',num2str(j)]}, ...
            'ResponseName', 'Y','FitPosterior',1, ...
            'ClassNames', labelconditions);
        Ypre(n) = predict(Mdl,Xtest);
        fprintf('*')
    end
    fprintf('\n')
    validationAccuracyOK(ns)=sum(Ypre'==Y)/numel(Y);
end
fprintf('\n')
AccuracySVM=mean(validationAccuracyOK);
% [mean(validationAccuracyOK),std(validationAccuracyOK)]
% PARAMETERS ############################################################
Mdl = fitcecoc(Xpca(:,[ONEpca,TWOpca]),Y, 'Learners', template, 'Coding',...
        'binarycomplete', 'PredictorNames', {['PCA_',num2str(ONEpca)],['PCA_',num2str(TWOpca)]}, ...
        'ResponseName', 'Y','FitPosterior',1, ...
        'ClassNames', labelconditions);
% Estimation of the labels:    
fprintf('>>Predicting: ')
[Yhat,~,~,pPosterior]=resubPredict(Mdl);
fprintf('done.\n')
%% Plots ROC & Confusion Matrix ##########################################
%Set binary targets
Targetss=zeros(Nobser,numel(unique(Y)));
Outputss=zeros(Nobser,numel(unique(Y)));
Yconditions=labelconditions;
for i=1:Nobser
    Targetss(i,Yconditions==Y(i))=1;
    Outputss(i,Yconditions==Yhat(i))=1;
end
% PLOT ROC curve
figure; plotroc(Targetss',Outputss');
ROCaxis=gca;
auxl=1;
for l=1:numel(condition_labels)
    ROCaxis.Children(auxl).Color=CM(ColorIndx(numel(condition_labels)-l+1),:);
    ROCaxis.Children(auxl).DisplayName=condition_labels{numel(condition_labels)-l+1};
    auxl=auxl+2;
end
% PLOT CONFUSION MATRIX
figure; plotconfusion(Targetss',Outputss')
CMaxis=gca;
for l=1:numel(condition_labels)
    CMaxis.YTickLabel{l}=condition_labels{l};
    CMaxis.XTickLabel{l}=condition_labels{l};
end
fprintf('\n%s\n',repmat('*',40,1))
disp(['LOOCV Acc SVM kernel ',kernels2test{kernelindx},...
    ' pC=',num2str(round(1000*AccuracySVM)/10),...
    '% +/-',num2str(num2str(round(1000*std(validationAccuracyOK))/10))]);
fprintf('\n%s\n',repmat('*',40,1))
%% POSTERIOR PROBABILITIES
deltaPlus=0.15;
MaxFactorPlus=ones(1,2);
MinFactorPlus=MaxFactorPlus;
maxSigns=sign(max(Xpca(:,[ONEpca,TWOpca])));
minSigns=sign(min(Xpca(:,[ONEpca,TWOpca])));
MaxFactorPlus(maxSigns>0)=1+deltaPlus;
MaxFactorPlus(maxSigns<0)=1-deltaPlus;
MinFactorPlus(minSigns>0)=1-deltaPlus;
MinFactorPlus(minSigns<0)=1+deltaPlus;
xMax = max(Xpca(:,[ONEpca,TWOpca])).*MaxFactorPlus;
xMin = min(Xpca(:,[ONEpca,TWOpca])).*MinFactorPlus;

x1Pts = linspace(xMin(1),xMax(1));
x2Pts = linspace(xMin(2),xMax(2));
[x1Grid,x2Grid] = meshgrid(x1Pts,x2Pts);
fprintf('>>Predicting Regions: ... ')
[~,~,~,PosteriorRegion] = predict(Mdl,[x1Grid(:),x2Grid(:)]);
fprintf('done.\n')
% Contour Probabilities
fprintf('>>Plotting eegions: ... ')
figure;
contourf(x1Grid,x2Grid,...
        reshape(max(PosteriorRegion,[],2),size(x1Grid,1),size(x1Grid,2)));
% drawnow
Contorno=gca;
% COLORMAP of the PROBAILITIES
% A=histcounts(PosteriorRegion);
CMprobs=cbrewer('seq','GnBu',20);
CMprobs=CMprobs(10:-1:1,:);
colormap(Contorno,CMprobs);
h = colorbar;
h.YLabel.String = 'Maximum posterior';
h.YLabel.FontSize = 8;
hold on
% SCATTER PLOT  WITH DATA POINTS
gh = gscatter(Xpca(:,ONEpca),Xpca(:,TWOpca),Y,'krb','ooo',12);
for n=1:numel(ColorIndx)
    gh(n).MarkerEdgeColor=[0,0,0];
    gh(n).MarkerFaceColor=CM(ColorIndx(n),:);
    gh(n).LineWidth=1;
end

title(['SCV kernel ',kernels2test{kernelindx},' LOOCV ACC=',num2str(round(1000*AccuracySVM)/10),'% ']);
xlabel([num2str(ONEpca),'th PCA component']);
ylabel([num2str(TWOpca),'th PCA component']);
axis([xMin(1),xMax(1),xMin(2),xMax(2)])
legend(gh,'Location','NorthWest')
hold off
fprintf(' done.\n')
%% FEATUREs CONTRIBUTION BIPLOT *******************************************
% A biplot allows you to visualize the magnitude and sign of each variable's
% contribution to the first two or three principal components, and how each 
% observation is represented in terms of those components.
figure;
hbi=biplot(coefforth(:,[ONEpca,TWOpca]),'scores',score(:,[ONEpca,TWOpca]),...
    'varlabels',NamesFeatures,'Marker','o','MarkerEdgeColor','k','MarkerSize',7);
% Vector Lines
for i=1:size(X,2)
    hbi(i).LineWidth=2;
end
% Vector Marker
for i=size(X,2)+1:2*size(X,2)
    hbi(i).MarkerFaceColor='b';
end
% Text Vector 
for i=2*size(X,2)+1:3*size(X,2)
    hbi(i).FontName='Arial';
    hbi(i).Interpreter='none';
end   
% Observations 
for i=3*size(X,2)+1:numel(hbi)-1
    hbi(i).MarkerFaceColor=CM(ColorIndx(labelconditions==Y(i-3*size(X,2))),:);
    hbi(i).MarkerSize=10;
end   

title('Feature contribution to each copmponents')
xlabel([num2str(ONEpca),' th PC']);
ylabel([num2str(TWOpca),' th PC']);
% axis([-.26 0.6 -.51 .51]);

%%  One-versus-all (OVA) Models ################

numClasses = numel(condition_labels);
inds = cell(numClasses,1);      % Preallocation
SVMModel = cell(numClasses,1);

% rng(1); % For reproducibility
fprintf('>>Making Binary Models: ')
for j = 1:numClasses
    Ybin= false(Nobser,1);
    Ybin(Y==condition_labels{j})=true;  % OVA classification
    SVMModel{j} = fitcsvm(Xpca(:,[ONEpca,TWOpca]),Ybin,'ClassNames',[false true],...
        'Standardize',true,'KernelFunction',kernels2test{kernelindx},...
        'KernelScale','auto');
    %         'KernelScale','auto');
    fprintf('%i,',j)
end
fprintf('\n')

fprintf('>>Calssifying Model: ')
for j = 1:numClasses
    %SVMModel{j} = fitPosterior(SVMModel{j});
    SVMModel{j} = fitSVMPosterior(SVMModel{j});
    fprintf('%i,',j)
end
fprintf('\n')

fprintf('>>Making Posterior Probability Regions: ')
% Same as Multiclass case
% d = 0.02;
% [x1Grid,x2Grid] = meshgrid(1.1*min(Xpca(:,ONEpca)):d:1.1*max(Xpca(:,ONEpca)),...
%     1.1*min(Xpca(:,TWOpca)):d:1.1*max(Xpca(:,TWOpca)));
xGrid = [x1Grid(:),x2Grid(:)];

posterior = cell(3,1); 
for j = 1:numClasses
    [~,posterior{j}] = predict(SVMModel{j},xGrid);
    MinProbs(j)=min(posterior{j}(:,2));
    MaxProbs(j)=max(posterior{j}(:,2));
    fprintf('%i,',j)
end
RangeProbs=[min(MinProbs),max(MaxProbs)];
fprintf('\n')

%% Figure OVAs ************************************************************
fprintf('>>Making Plot: ...\n')
figure
NcolsFig=nextpow2(numClasses+1);
haxis = zeros(numClasses + 1,1); % Preallocation for graphics handles
for j = 1:numClasses
    fprintf('>>Surface of Model %i: ...',j)
    Ax=subplot(2,NcolsFig,j);
    contourf(x1Grid,x2Grid,reshape(posterior{j}(:,2),size(x1Grid,1),size(x1Grid,2)));
    colormap(CMprobs);
    caxis(RangeProbs);
    hold on
    % h(1:numClasses) = gscatter(Xpca(:,ONEpca),Xpca(:,TWOpca),Y,CM(ColorIndx,:),'ooo',10);
    h = gscatter(Xpca(:,ONEpca),Xpca(:,TWOpca),Y,CM(ColorIndx,:),'ooo',8);
    for n=1:numel(ColorIndx)
        h(n).MarkerEdgeColor=[0,0,0];
        h(n).MarkerFaceColor=CM(ColorIndx(n),:);
        h(n).LineWidth=0.5;
    end
    title(sprintf('Posteriors for %s Class',condition_labels{j}));
    xlabel([num2str(ONEpca),' th PCA component']);
    ylabel([num2str(TWOpca),' th PCA component']);
    legend off
    axis tight
    hold off
    fprintf('done.\n')
end

% hcm=colorbar;
% hcm.YLabel.String = 'Maximum posterior';
% hcm.YLabel.FontSize = 8;
haxis(numClasses + 1) = colorbar('Location','EastOutside',...
    'Position',[[0.8,0.1,0.05,0.4]]);

set(get(haxis(numClasses + 1),'YLabel'),'String','Posterior','FontSize',9);
legend(h(1:numClasses),'Location',[0.6,0.2,0.1,0.1]);
fprintf('\n\n>>[PCA->SVM to Neuronal Features]: Ready.\n\n')
%% Save Results:
% Data:     X,Y [Previously done]
% Models:
%   SVMModel  OVA models
%   Mdl       Multiclass
checkname=1;
while checkname==1
    % Get Directory
    fprintf('\n\n [ To add models, please select the .mat file ] \n\n')
    [FileName,PathName] = uigetfile('*.mat',[' Pick the Data table File to ADD models: ',FileMat],...
        'MultiSelect', 'off',FolderName);
    % dotindex=find(FileName=='.');
    if strcmp(FileName,FileMat)
        checkname=0;
        % SAVE DATA
        save([PathName,FileName],'SVMModel','Mdl','-append');
        disp([FileMat,'   -> MODELS ADDED '])
    elseif FileName==0
        checkname=0;
        disp('No model added')
    else
        disp('Not the same File!')
        disp('Try again!')
    end
end  
%% END OF THE WORLD *******************************************************