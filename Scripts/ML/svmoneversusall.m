%%  One-versus-all (OVA) Models ################
% Build model using : one-leave-out
% Input
%   Xpca,               Matrix of PCs
%   ONEpca,             PC i
%   TWOpca,             PC j
%   Y,                  Vector of Labels of K unqiue labels
%   BestKernel,         Kernel for SVM
%   x1Grid,             Values to evaluate of PC i
%   x2Grid,             Values to evaluate of PC i
% Output
%   posterior           Grid of probabilities
%   Yallbin             Cell of predicted for each {Class}(Nobs,n-model)
%   RangeProbs          Range of Values for Colormap probabilities
% 
function [ProbSups,Yallbin,RangeProbs,Ybestova]=svmoneversusall(Xpca,ONEpca,TWOpca,Y,BestKernel,Nsim,x1Grid,x2Grid)
% Setup 
labelconditions=unique(Y,'legacy');
Nobser=size(Xpca,1);
numClasses = numel(labelconditions);
xGrid = [x1Grid(:),x2Grid(:)]; % Grid Values PC i-th vs PC j-th
% Models
SVMModel = cell(numClasses,1);
% Binary Outputs
Yallbin= cell(numClasses,1);
% Grid of probabilities
posterior = cell(3,1);
tablecount=tabulate(Y);
Noun=zeros(numClasses,1);
% Intialize
Ybestova=categorical(zeros(Nobser,Nsim));
for j = 1:numClasses
    Noun(j)=tablecount{j,2};
    ProbSups{j}=zeros(size(xGrid,1),1);
    Yallbin{j}=zeros(Nobser,Nsim);
end
kfolder=min(Noun);
fprintf('>>Making Binary OVA Models: \n')
for n=1:Nsim
    % SVM model ofr each class:
    fprintf('>>Model SVM one-versus-all %i / %i: ',n,Nsim);
    Yova=categorical(zeros(Nobser,1)); % categorical predictions of K models
    % Class SVM models
    ProbsY=zeros(Nobser,numClasses);
    % Ybinclasses=zeros(Nobser,numClasses);
    for j = 1:numClasses
        Ybin= false(Nobser,1); % init biniry target s
        Ybin(Y==labelconditions(j))=true;  % OVA classification
        SVMModel{j} = fitcsvm(Xpca(:,[ONEpca,TWOpca]),Ybin,'ClassNames',[false true],...
            'Standardize',true,'KernelFunction',BestKernel,...
            'KernelScale','auto');
        SVMModel{j} = fitSVMPosterior(SVMModel{j},'KFold',kfolder);
        % SVMModel{j} = fitSVMPosterior(SVMModel{j},'Holdout',0.2);
        % SVMModel{j} = fitSVMPosterior(SVMModel{j},'Leaveout','on');
        [Yhatbin,pY] = predict(SVMModel{j},Xpca(:,[ONEpca,TWOpca]));
        % py: probs of negative and positive class: p([-,+])
        ProbsY(:,j)=pY(:,2); % posterior prob of j model
        Yova(Yhatbin>0,1)=labelconditions(j); % make it categorical
        % Binary output for j-th class and n-th model
        Yallbin{j}(:,n)=Yhatbin;
        % Ybinclasses(:,j)=Yhatbin; % Multiclass output
        % Grid values evaluation
        [~,posterior{j}] = predict(SVMModel{j},xGrid);
        ProbSups{j}=ProbSups{j}+posterior{j}(:,2);
        fprintf('%i,',j)
    end
%     % Look for cases of confused labels if any [not really necessary]
%     Yaux=sum(Ybinclasses,2);
%     LabConf=find(Yaux>1); % Confused Rows
%     if isempty(LabConf)
%         fprintf('>No Confused Labels found\n');
%     else
%         fprintf('>Confused Labels:\n');
%         for k=1:numel(LabConf)
%             labsconf=find(Ybinclasses(LabConf(k),:)==true);
%             fprintf('>');
%             for i=1:numel(labsconf)
%                 fprintf('%s (%3.2f)',char(labelconditions(labsconf(i))),...
%                     ProbsY(LabConf(k),labsconf(i)))
%                 if i<numel(labsconf) fprintf(' & '); end
%             end
%             [~,oklab]=max(ProbsY(LabConf(k),:));
%             rejectlabs=setdiff(1:numClasses,oklab);
%             Ybinclasses(LabConf(k),rejectlabs)=false;
%             fprintf('\n')
%         end
%         fprintf('\n')
%     end
    % After Building the K models, choose more likely label
    [~,IndexClass]=max(ProbsY'); 
    for k=1:Nobser
        Ybestova(k,n)=labelconditions(IndexClass(k));
    end
end
% Output details:
for j = 1:numClasses
    ProbSups{j}=ProbSups{j}/Nsim;
    MinProbs(j)=min(ProbSups{j});
    MaxProbs(j)=max(ProbSups{j});
end
RangeProbs=[min(MinProbs),max(MaxProbs)];
fprintf('\n')