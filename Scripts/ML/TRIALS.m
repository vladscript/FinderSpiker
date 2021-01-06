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
        % LEAVE ONE OUT
        
        fprintf('\n>>Model %i/%i,kernel %s, PCs: %i vs %i -> \n',c,Ncombo,kernels2test{k},i,j);
        for ns=1:Nsim
            for n=1:Nobser
                Xtrain=Xdata(setdiff(1:Nobser,n),:);
                Ytrain=Y(setdiff(1:Nobser,n),:);
                Xtest=Xdata(n,:);
                Ytest=Y(n,:);
                Mdl = fitcecoc(Xtrain,Ytrain, 'Learners', template, 'Coding',...
                    'binarycomplete', 'PredictorNames', {['PCA_',num2str(i)],['PCA_',num2str(j)]}, ...
                    'ResponseName', 'Y','FitPosterior',1, ...
                    'ClassNames', labelconditions);
                Yhatsvm(n) = predict(Mdl,Xtest);
                fprintf('*')
            end
            fprintf('\n')
            validationAccuracyOK(ns)=sum(Yhatsvm'==Y)/numel(Y);
        end
        LOOCVs{c,k}=validationAccuracyOK;
        ALLPCAs(c,2+k)=mean(validationAccuracyOK);
        fprintf('\n>>Mean LOOCV accuracy: %3.2f\n',mean(validationAccuracyOK));
%         % SVM
%         Mdl = fitcecoc(Xpca(:,[i,j]),Y, 'Learners', template, 'Coding',...
%         'binarycomplete', 'PredictorNames', {['PCA_',num2str(i)],['PCA_',num2str(j)]}, ...
%         'ResponseName', 'Y','FitPosterior',1, ...
%         'ClassNames', labelconditions);
%         %[Yhat,~,~,pPosterior]=resubPredict(Mdl);
%         % Evlauate MODEL: Cross Validation Methods
%         % partitionedModel = crossval(Mdl, 'KFold', 5); % not good
%         partitionedModel = crossval(Mdl,'leaveout','on'); % low number of observations
%         validationAccuracyOK = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
%         ALLPCAs(c,2+k)=validationAccuracyOK;
%         [Yhat,~,~,~]=resubPredict(Mdl);
%         validationAccuracy=sum(Y==Yhat)/numel(Y);
        % Save results
        % PCAandACC=[PCAandACC;i,j,validationAccuracy];
        
    end
end
fprintf(repmat('*',10,1)); fprintf('[complete]\n');
% GET BEST Kernel multiclassificator
if size(ALLPCAs(:,3:end),1)>1
    MaXeachKernel=max(ALLPCAs(:,3:end));
else
    MaXeachKernel=ALLPCAs(:,3:end);
end
[AccuracySVM,kernelindx]=max(MaXeachKernel);

if kernelindx<3
    template = templateSVM('KernelFunction',kernels2test{kernelindx}, 'PolynomialOrder', [],...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);    
%     'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);    
else
    template = templateSVM('KernelFunction',kernels2test{kernelindx}, 'PolynomialOrder', 2,...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);    
%     'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);
end
% BEST Model and PCAs *************************
[~,nModel]=max(ALLPCAs(:,2+kernelindx));
columns2boxplot(LOOCVs{nModel,1}',LOOCVs{nModel,2}',LOOCVs{nModel,3}',{'Lineal','Gaussian','Polynomial'})
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

%% #############################################################

MdlSVM = fitcecoc(Xpca(:,[ONEpca,TWOpca]),Y, 'Learners', template, 'Coding',...
        'binarycomplete', 'PredictorNames', {['PCA_',num2str(ONEpca)],['PCA_',num2str(TWOpca)]}, ...
        'ResponseName', 'Y','FitPosterior',1, ...
        'ClassNames', labelconditions);
    
Mdl.BinaryLearners{1}.KernelParameters.Scale
Mdl.BinaryLearners{2}.KernelParameters.Scale
Mdl.BinaryLearners{3}.KernelParameters.Scale
Mdl.LearnerWeights

MdlSVM.BinaryLearners{1}.KernelParameters.Scale


for n=1:Nobser
                Xtrain=Xdata(setdiff(1:Nobser,n),:);
                Ytrain=Y(setdiff(1:Nobser,n),:);
                Xtest=Xdata(n,:);
                Ytest=Y(n,:);
                Mdl = fitcecoc(Xtrain,Ytrain, 'Learners', template, 'Coding',...
                    'binarycomplete', 'PredictorNames', {['PCA_',num2str(i)],['PCA_',num2str(j)]}, ...
                    'ResponseName', 'Y','FitPosterior',1, ...
                    'ClassNames', labelconditions);
                Yhatsvm(n) = predict(Mdl,Xtest);
                fprintf('*')
            end