%% PCA & SVM Classification of
% Feature Matrix:   X
% Label vector:     Y
% PCA: variance weighted
% SVM:
%   Multiclass:     fitcecoc
%   One-Versus-All: fitcsvm
% Read Table(s) and Select Features from Datasets:
%  - Raster_Activity_Dataset_XXXX.csv
%  - General_Ensembles_Dataset_YYYY.CSV
%  - Network_Features_Dataset_ZZZZ.csv
%  - Or CSV tables with Features/ Columns and ID columns
% IMPORTANT:
% It assumes IDs are the same among tables (for sorting purposes)
% 
%% Setup:
clc; clear; % close all;
% Import_FinderSpiker;
Nsim=20; % Repetitions for SVM
% Probaility colormap:
CMprobs=cbrewer('seq','GnBu',20);
CMprobs=CMprobs(10:-1:1,:);
% 
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
        NaNExps=find(isnan(sum(X,2)));
        OkExps=find(~isnan(sum(X,2)));
        X=X(OkExps,:);
        Y=Y(OkExps);
        NaNExps=find(isnan(sum(X,2)));
        
end

%% Color for Labels 
[Nobser,Nfeatures]=size(X);
labelconditions=unique(Y,'legacy');

for n=1:numel(labelconditions)
    condition_labels{n}=char(labelconditions(n));
end

[CM,ColorIndx]=Color_Selector(condition_labels);
%% Plot Regions (?)
% It all starts with a choice ...
choice = questdlg('Plot classification regions?', ...
	'Probability of PC on SVM', ...
	'Plot Regions','Do not plot','Cancel','Plot Regions');
%% PCA 
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
addpc=false;
if Npcs==1 Npcs=2; addpc=true; end
% Only the first 95% of the cumulative distribution is displayed
figure; 
pareto(explained); hold on
xlabel('Principal Component')
ylabel('Variance Explained (%)')
% plot([Npcs,Npcs],[0,100],'r');

Xpca=score(:,1:Npcs);

%% Principal Components & Kernel Selection ###############################
% Features:     Xpca,   dims:   Nobservations x Npcs
% Labels:       Y,      dims:   Nobservations x 1
fprintf('\n>> For every kernel and pair of PCs get best accuracy: ')
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
        fprintf('>>Model %i/%i,kernel %s, PCs: %i vs %i -> \n',c,Ncombo,kernels2test{k},i,j);
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
        fprintf('>>Mean Accuracy: %3.2f\n',mean(validationAccuracyOK));        
    end
end
% fprintf(repmat('*',10,1)); fprintf('[complete]\n');
% GET BEST Kernel multiclassificator
if size(ALLPCAs(:,3:end),1)>1
    MaXeachKernel=max(ALLPCAs(:,3:end));
else
    MaXeachKernel=ALLPCAs(:,3:end);
end
[~,kernelindx]=max(MaXeachKernel);
BestKernel=kernels2test{kernelindx};
if kernelindx<3
    template = templateSVM('KernelFunction',BestKernel, 'PolynomialOrder', [],...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);    
%     'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);    
else
    template = templateSVM('KernelFunction',BestKernel, 'PolynomialOrder', 2,...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);    
%     'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);
end
% BEST PCAs *************************
[~,nModel]=max(ALLPCAs(:,2+kernelindx));
% PCs Matrix Average
matPCs=zeros(Npcs);
for n=1:Ncombo
    i=ALLPCAs(n,1);
    j=ALLPCAs(n,2);
    matPCs(i,j)=ALLPCAs(n,2+kernelindx);
    matPCs(j,i)=matPCs(i,j);
end
columns2boxplot(LOOCVs{nModel,1}',LOOCVs{nModel,2}',LOOCVs{nModel,3}',{'Lineal','Gaussian','Polynomial'})
ONEpca=ALLPCAs(nModel(1),1);
TWOpca=ALLPCAs(nModel(1),2);
%% Best Classifying PC Pair contribution BIPLOT 
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
%% Leave-One-Out Cross Validation OF THE MODEL 
% clear Yhatsvm
fprintf('Model| Leave-One-Out Cross Validations \n')
Ysvmall=categorical(zeros(Nobser,Nsim));
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
    Ysvmall(:,ns)=Ypre';
    fprintf('\n')
    validationAccuracyOK(ns)=sum(Ypre'==Y)/numel(Y);
end
fprintf('\n')
AccuracySVM=mean(validationAccuracyOK);
[ROCloocv,METRICSloocv]=rockandconfusion(Ysvmall,Y);

%% Plot Classificaion Regions #############################################
switch choice
    case 'Plot Regions'
        %% Multiclass Regions binary complete
        fprintf('\n>>Making PCs grid: ')
        deltaPlus=0.15; % Resolution of grid/contour
        [x1Grid,x2Grid]=makepcagrid(Xpca,ONEpca,TWOpca,deltaPlus);
        % [x1Grid,x2Grid]=makepcagrid(Xpca,ONEpca,TWOpca,deltaPlus,200); %ALSO
        fprintf('done.\n')
        % MULTICLASS MODELS
        [PosteriorRegion,Ymany]=svmmulticlass(Xpca,ONEpca,TWOpca,Y,template,labelconditions,Nsim,x1Grid,x2Grid);
        % Get securely predicted
        [Yclass,ClassRateMC]=getrightclass(Ymany,Y);
        % PLOT Regions for Multiclass
        plot_multiclass(x1Grid,x2Grid,PosteriorRegion,CMprobs,Xpca,ONEpca,...
            TWOpca,Y,CM,ColorIndx,BestKernel,AccuracySVM,Yclass);
        % ROCs & F1 Scores
        [ROC,METRICS]=rockandconfusion(Ymany,Y);
        % Get MissClassified
        %% One-versus-all Regions
        [posterior,Yallbin,RangeProbs,Ybestova]=svmoneversusall(Xpca,ONEpca,TWOpca,Y,BestKernel,Nsim,x1Grid,x2Grid);
        % Get right classes:
        [Yconfident,ClassRatebin]=getrighclassbinary(Y,Yallbin);
        % PLOT Regions for OVA
        plot_ova(x1Grid,x2Grid,posterior,CMprobs,Xpca,ONEpca,labelconditions,...
            TWOpca,Y,CM,ColorIndx,RangeProbs,Yconfident);
        % MAKE THEM CATEGORICAL EACH CLASS AND GET METRICS
         for k=1:numel(labelconditions)
             Currentlab=labelconditions(k);
             % Predicted
             PredBin=Yallbin{k};
             YmanyOVA=categorical(false(Nobser,Nsim));
             for n=1:Nsim
                 YmanyOVA(PredBin(:,n)>0,n)=Currentlab;
             end
             % Gorund - truth: class-false
             Yovabin=categorical(false(Nobser,1));
             Yovabin(Y==labelconditions(k))=labelconditions(k); 
             % ROCs & F1 Scores
            [ROCova{k},METRICSova{k}]=rockandconfusion(YmanyOVA,Yovabin);
         end

    case 'Do not plot'
        fprintf('\n>>End.\n')
    case 'Cancel'
        fprintf('\n>>End.\n')
end


%% Resume Data
fprintf('****************************************\n');
fprintf('PCA and SVM classification of:\n>> %i features, labels:\n', Nfeatures)
tabulate(Y);
fprintf('  PCs with >95%% of data variance: ');
if addpc 
    fprintf('%i\n',Npcs-1);
else
    fprintf('%i\n',Npcs);
end
fprintf('Best Accuracy with %s kernel and %i and %i PCs',...
    BestKernel,ONEpca,TWOpca);
fprintf('\nMean LOOCV Accuracy (SD): %3.2f%%(%3.2f)\n',100*mean(validationAccuracyOK),std(100*validationAccuracyOK))
% Make table of mean  and standard deviation of all metrics in LOOCV
% For Multiclass binary complete vs one versus all
FinalTableMC=metrotable(METRICS,labelconditions);
for c=1:numel(labelconditions)
    METRICStrue{c}=METRICSova{c}{2};
end
FinalTableOVA=metrotable(METRICStrue,labelconditions);
fprintf('****************************************');

% %% Show Classification Metrics
% Y=DataXY.Label;
% labelconditions=unique(Y,'legacy');
% 
% FinalTableMC=metrotable(METRICS,labelconditions);
% for c=1:numel(labelconditions)
%     METRICStrue{c}=METRICSova{c}{2};
% end
% FinalTableOVA=metrotable(METRICStrue,labelconditions);
% disp('OK');

%% Save Results: #########################################################
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
        save([PathName,FileName],'ALLPCAs','Xpca','Npcs','ONEpca',...
            'BestKernel','TWOpca','METRICSloocv','METRICS','METRICSova',...
            '-append');
        disp([FileMat,'   -> METRICS SAVED '])
    elseif FileName==0
        checkname=0;
        disp('No model added')
    else
        disp('Not the same File!')
        disp('Try again!')
    end
end  
%% END OF THE WORLD