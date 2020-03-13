%% Get Data ##############################################################
% Read Table and Select Features

% X & Y From Feature Tables [in progress--]

[Nobser,Nfeatures]=size(X);
% Interface to merge conditions
% Merge Same Condition/Different Experiments
Y(Y=='DyskinesiaA')='Dyskinesia';
Y(Y=='DyskinesiaC')='Dyskinesia';


Y = removecats(Y); % removes empty categories
labelconditions=unique(Y);
condition_labels={};
for i=1:numel(labelconditions)
    fprintf('+')
    condition_labels{i}=char(labelconditions(i));
end
[CM,ColorIndx]=Color_Selector(condition_labels);

%% PCA ###################################################################
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
% Display 2 PCA components scatter plot:
figure
gscatter(score(:,1),score(:,4),Y,CM(ColorIndx,:),'dso',10);
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
% Rows of score correspond to observations, 
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
% Test 3 function kernels: linear, parabolic and gaussian
kernels2test={'linear';'gaussian';'polynomial'};
Combo=combnk(1:Npcs,2);
Ncombo=size(Combo,1);
ALLPCAs=[Combo,zeros(Ncombo,3)];
for k=1:numel(kernels2test)
    % generate SVM template
    fprintf('\n>>Making SVM Template with kernel: %s ',kernels2test{k})
    if k<3
        template = templateSVM('KernelFunction',kernels2test{k}, 'PolynomialOrder', [],...
        'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);    
    else
        template = templateSVM('KernelFunction',kernels2test{k}, 'PolynomialOrder', 2,...
        'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',0);
        fprintf(' of order : 2')
    end
    fprintf('\n')
    % SVM Classificator using 2 PCAs: all combinations
    for c=1:Ncombo
        i=ALLPCAs(c,1); j=ALLPCAs(c,2);
        % SVM
        Mdl = fitcecoc(Xpca(:,[i,j]),Y, 'Learners', template, 'Coding',...
        'binarycomplete', 'PredictorNames', {['PCA_',num2str(i)],['PCA_',num2str(j)]}, ...
        'ResponseName', 'Y','FitPosterior',1, ...
        'ClassNames', labelconditions);
        %[Yhat,~,~,pPosterior]=resubPredict(Mdl);
        % Evlauate MODEL: Cross Validation Methods
        % partitionedModel = crossval(Mdl, 'KFold', 5); % not good
        partitionedModel = crossval(Mdl,'leaveout','on'); % low number of observations
        validationAccuracyOK = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
        ALLPCAs(c,2+k)=validationAccuracyOK;
%         [Yhat,~,~,~]=resubPredict(Mdl);
%         validationAccuracy=sum(Y==Yhat)/numel(Y);
        % Save results
        % PCAandACC=[PCAandACC;i,j,validationAccuracy];
        fprintf(repmat('*',30,1));
        fprintf('\n>>Model %i/%i,kernel %s, PCs: %i vs %i -> %3.2f\n',c,Ncombo,kernels2test{k},i,j,validationAccuracyOK);
        fprintf(repmat('*',20,1)); fprintf('[classifying]');
    end
end
fprintf(repmat('*',10,1)); fprintf('[complete]\n');
% GET BEST Kernel multiclassificator
MaXeachKernel=max(ALLPCAs(:,3:end));
[AccuracySVM,kernelindx]=max(MaXeachKernel);

if kernelindx<3
    template = templateSVM('KernelFunction',kernels2test{kernelindx}, 'PolynomialOrder', [],...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);    
else
    template = templateSVM('KernelFunction',kernels2test{kernelindx}, 'PolynomialOrder', 2,...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);
end
% BEST Model and PCAs
[~,nModel]=max(ALLPCAs(:,2+kernelindx));
ONEpca=ALLPCAs(nModel(1),1);
TWOpca=ALLPCAs(nModel(1),2);
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

%% POSTERIOR PROBABILITIES
xMax = max(Xpca(:,[ONEpca,TWOpca]));
xMin = min(Xpca(:,[ONEpca,TWOpca]));

x1Pts = linspace(xMin(1),xMax(1));
x2Pts = linspace(xMin(2),xMax(2));
[x1Grid,x2Grid] = meshgrid(x1Pts,x2Pts);
fprintf('>>Predicting Regions: ... ')
[~,~,~,PosteriorRegion] = predict(Mdl,[x1Grid(:),x2Grid(:)]);
fprintf('>>Done.\n')
figure;
contourf(x1Grid,x2Grid,...
        reshape(max(PosteriorRegion,[],2),size(x1Grid,1),size(x1Grid,2)));
h = colorbar;
h.YLabel.String = 'Maximum posterior';
h.YLabel.FontSize = 15;
hold on
gh = gscatter(Xpca(:,ONEpca),Xpca(:,TWOpca),Y,'krb','*xd',8);
for n=1:numel(ColorIndx)
    gh(n).Color=CM(ColorIndx(n),:);
    gh(n).LineWidth=2;
end

title(['SCV kernel ',kernels2test{kernelindx},' pC=',num2str(round(1000*AccuracySVM)/10),'% & Maximum Posterior']);
xlabel([num2str(ONEpca),'th PCA component']);
ylabel([num2str(TWOpca),'th PCA component']);
axis tight
legend(gh,'Location','NorthWest')
hold off
%% GSCATTERS per CONDITION one-versus-all (OVA) ################

numClasses = numel(condition_labels);
inds = cell(numClasses,1);      % Preallocation
SVMModel = cell(numClasses,1);

% rng(1); % For reproducibility
fprintf('>>Making Binary Models: ')
for j = 1:numClasses
    Ybin= false(Nobser,1);
    Ybin(Y==condition_labels{j})=true;  % OVA classification
    SVMModel{j} = fitcsvm(Xpca(:,[ONEpca,TWOpca]),Ybin,'ClassNames',[false true],...
        'Standardize',true,'KernelFunction','rbf','KernelScale','auto');
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

d = 0.02;
[x1Grid,x2Grid] = meshgrid(min(Xpca(:,ONEpca)):d:max(Xpca(:,ONEpca)),...
    min(Xpca(:,TWOpca)):d:max(Xpca(:,TWOpca)));
xGrid = [x1Grid(:),x2Grid(:)];

posterior = cell(3,1); 
for j = 1:numClasses
    [~,posterior{j}] = predict(SVMModel{j},xGrid);
    fprintf('%i,',j)
end
fprintf('\n')

fprintf('>>Making Plot: ...\n')
figure
h = zeros(numClasses + 1,1); % Preallocation for graphics handles
for j = 1:numClasses
    fprintf('>>Surface of Model %i: ...',j)
    subplot(2,2,j)
    contourf(x1Grid,x2Grid,reshape(posterior{j}(:,2),size(x1Grid,1),size(x1Grid,2)));
    hold on
    h(1:numClasses) = gscatter(Xpca(:,ONEpca),Xpca(:,TWOpca),Y,CM(ColorIndx,:),'*xd',10);
    title(sprintf('Posteriors for %s Class',condition_labels{j}));
    xlabel([num2str(ONEpca),' th PCA component']);
    ylabel([num2str(TWOpca),' th PCA component']);
    legend off
    axis tight
    hold off
    fprintf('done.\n')
end
h(numClasses + 1) = colorbar('Location','EastOutside',...
    'Position',[[0.8,0.1,0.05,0.4]]);
set(get(h(numClasses + 1),'YLabel'),'String','Posterior','FontSize',16);
legend(h(1:numClasses),'Location',[0.6,0.2,0.1,0.1]);
fprintf('\n\n>>[PCA->SVM to Neuronal Features]: Ready.\n\n')
%% Save Results:
% Data:     X,Y
% Models:   Finest Multiclass SVM and SVM models
% Plots:    PCA,MC,ROC,MultiClassRegions,OnevsAllRegions
%% END OF THE WORLD *******************************************************