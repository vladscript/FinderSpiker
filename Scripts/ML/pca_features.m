%% Get Data ##############################################################

% X=copy data from ALLDATA.xlsx or Feature Tables [in progress--]

% Merge Same Condition/Different Experiments

Y(Y=='DyskinesiaA')='Dyskinesia';
Y(Y=='DyskinesiaC')='Dyskinesia';

% maybe an interface to merge conditions ... (but how?)

Y = removecats(Y); % removes empty categories
labelconditions=unique(Y);
for i=1:numel(labelconditions)
    condition_labels{i}=char(labelconditions(i));
end
[CM,ColorIndx]=Color_Selector(condition_labels);

%% PCA ###################################################################
% function Xpca=pca_features(X,varLevel)
fprintf('>>PCA :')
[coefs,score,latent] = pca(X,'algorithm','als');
% Rows of score correspond to observations, 
% and columns to components.
VarExplained=cumsum(latent)./sum(latent);
% Npcs=find(VarExplained>=0.95,1);
Xpca=coefs*score'; % Features c observations
[Nfeatures,Nobser]=size(Xpca);

fprintf('done\n')

%% Multiclass Classification #############################################
% Test 3 function kernels: linear, parabolic and gaussian
kernels2test={'linear';'gaussian';'polynomial'};
ALLPCAs=[];
for k=2:3
    % generate SVM template
    fprintf('>>Making SVM Template with kernel: %s ',kernels2test{k})
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
    aux=0; PCAandACC=[];
    for i=1:Nfeatures-1
        for j=i+1:Nfeatures
            % SVM
            Mdl = fitcecoc(Xpca([i,j],:)',Y, 'Learners', template, 'Coding',...
            'onevsone', 'PredictorNames', {'Var1_1' 'Var1_9'}, ...
            'ResponseName', 'Y','FitPosterior',1, ...
            'ClassNames', labelconditions);
            %[Yhat,~,~,pPosterior]=resubPredict(Mdl);
            % Evlauate 
            partitionedModel = crossval(Mdl, 'KFold', 5);
            validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
            % Save results
            PCAandACC=[PCAandACC;i,j,validationAccuracy];
            fprintf(repmat('*',20,1));
            fprintf('\n>>Model %i,kernel %s, PCA components: %i vs %i\n',aux,kernels2test{k},i,j);
            fprintf(repmat('*',20,1));
            aux=aux+1;
        end
    end
    % Ssave Best accuracy for each kernel
    nModelKernel(k)=max(PCAandACC(:,3));
    ALLPCAs(:,:,k)=PCAandACC;
end
% GET BEST Kernel multiclassificator
[AccuracySVM,kernelindx]=max(nModelKernel);
PCAandACC=ALLPCAs(:,:,kernelindx);
if kernelindx<3
    template = templateSVM('KernelFunction',kernels2test{kernelindx}, 'PolynomialOrder', [],...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);    
else
    template = templateSVM('KernelFunction',kernels2test{kernelindx}, 'PolynomialOrder', 2,...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);
end
% BEST Model and PCAs
[~,nModel]=max(PCAandACC(:,3));
ONEpca=PCAandACC(nModel(1),1);
TWOpca=PCAandACC(nModel(1),2);
Mdl = fitcecoc(Xpca([ONEpca,TWOpca],:)',Y, 'Learners', template, 'Coding',...
        'onevsone', 'PredictorNames', {'Var1_1' 'Var1_9'}, ...
        'ResponseName', 'Y','FitPosterior',1, ...
        'ClassNames', labelconditions);
% Estimation of the labels:    
[Yhat,~,~,pPosterior]=resubPredict(Mdl);

%% PLots For the Best;(ROC COnfision MAtrix, etc)
%Set binary targets
Targetss=zeros(Nobser,numel(unique(Y)));
Outputss=zeros(Nobser,numel(unique(Y)));
Yconditions=unique(Y,'legacy');
for i=1:Nobser
    Targetss(i,Yconditions==Y(i))=1;
    Outputss(i,Yconditions==Yhat(i))=1;
end
% PLOT ROC curve
figure; plotroc(Targetss',Outputss');
ROCaxis=gca;
auxl=1;
for l=1:numel(classNames)
    ROCaxis.Children(auxl).Color=CM(l,:);
    ROCaxis.Children(auxl).DisplayName=classNames{l};
    auxl=auxl+2;
end
% PLOT CONFUSION MATRIX
figure; plotconfusion(Targetss',Outputss')
CMaxis=gca;
for l=1:numel(classNames)
    CMaxis.YTickLabel{l}=classNames{l};
    CMaxis.XTickLabel{l}=classNames{l};
end

%% POSTERIOR PROBABILITIES
xMax = max(Xpca([ONEpca,TWOpca],:)');
xMin = min(Xpca([ONEpca,TWOpca],:)');

x1Pts = linspace(xMin(1),xMax(1));
x2Pts = linspace(xMin(2),xMax(2));
[x1Grid,x2Grid] = meshgrid(x1Pts,x2Pts);

[~,~,~,PosteriorRegion] = predict(Mdl,[x1Grid(:),x2Grid(:)]);

figure;
contourf(x1Grid,x2Grid,...
        reshape(max(PosteriorRegion,[],2),size(x1Grid,1),size(x1Grid,2)));
h = colorbar;
h.YLabel.String = 'Maximum posterior';
h.YLabel.FontSize = 15;
hold on
gh = gscatter(Xpca(ONEpca,:),Xpca(TWOpca,:)',Y,'krb','*xd',8);
gh(1).Color=CM(2,:);
gh(1).LineWidth=2;
gh(2).Color=CM(3,:);
gh(2).LineWidth=2;
gh(3).Color=CM(1,:);
gh(3).LineWidth=2;

title(['SCV kernel ',kernels2test{kernelindx},' pC=',num2str(round(1000*AccuracySVM)/10),'% & Maximum Posterior']);
xlabel([num2str(ONEpca),'th PCA component']);
ylabel([num2str(TWOpca),'th PCA component']);
axis tight
legend(gh,'Location','NorthWest')
hold off
%% GSCATTERS per CONDITION one-versus-all (OVA) ################

% figure
% gscatter(Xpca(ONEpca,:),Xpca(TWOpca,:),Y,CM([2,3,1],:),'dso',10);
% title('{\bf Best calssifyng PCA components}');
% xlabel([num2str(ONEpca),' th PCA component']);
% ylabel([num2str(TWOpca),' th PCA component']);
% legend('Location','Northwest'); 
% axis tight

classNames = {'Amantadine';'Clozapine';'Dyskinesia'};
numClasses = size(classNames,1);
inds = cell(3,1); % Preallocation
SVMModel = cell(3,1);

% rng(1); % For reproducibility
for j = 1:numClasses
    Ybin= false(Nobser,1);
    Ybin(Y==classNames{j})=true;  % OVA classification
    SVMModel{j} = fitcsvm(Xpca([ONEpca,TWOpca],:)',Ybin,'ClassNames',[false true],...
        'Standardize',true,'KernelFunction','rbf','KernelScale','auto');
end

for j = 1:numClasses
    SVMModel{j} = fitPosterior(SVMModel{j});
end

d = 0.002;
[x1Grid,x2Grid] = meshgrid(min(Xpca(ONEpca,:)):d:max(Xpca(ONEpca,:)),...
    min(Xpca(TWOpca,:)):d:max(Xpca(TWOpca,:)));
xGrid = [x1Grid(:),x2Grid(:)];

posterior = cell(3,1); 
for j = 1:numClasses
    [~,posterior{j}] = predict(SVMModel{j},xGrid);
end

figure
h = zeros(numClasses + 1,1); % Preallocation for graphics handles
for j = 1:numClasses
    subplot(2,2,j)
    contourf(x1Grid,x2Grid,reshape(posterior{j}(:,2),size(x1Grid,1),size(x1Grid,2)));
    hold on
    h(1:numClasses) = gscatter(Xpca(ONEpca,:),Xpca(TWOpca,:),Y,CM([2,3,1],:),'dso',10);
    title(sprintf('Posteriors for %s Class',classNames{j}));
    xlabel([num2str(ONEpca),' th PCA component']);
    ylabel([num2str(TWOpca),' th PCA component']);
    legend off
    axis tight
    hold off
end
h(numClasses + 1) = colorbar('Location','EastOutside',...
    'Position',[[0.8,0.1,0.05,0.4]]);
set(get(h(numClasses + 1),'YLabel'),'String','Posterior','FontSize',16);
legend(h(1:numClasses),'Location',[0.6,0.2,0.1,0.1]);

%% END OF THE WORLD

% %% RAINCLOUDS of THE BEST PCAs
% pcaAMAN=Xpca(ONEpca,Y=='Amantadine');
% pcaCLZ=Xpca(ONEpca,Y=='Clozapine');
% pcaDYSK=Xpca(ONEpca,Y=='Dyskinesia');
% 
% figure
% raincloud_plot(pcaAMAN,'color',CM(ColorIndx(2),:),'box_on',1,...
%     'alphaval',5,'box_dodge',1,...
%     'box_dodge_amount',0.4, 'dot_dodge_amount', 0.4,...
%     'box_col_match',0,'box_dodge_amount',.25, 'dot_dodge_amount', 0.25,...
%     'line_width',3,'lwr_bnd',2);
% 
% raincloud_plot(pcaCLZ,'color',CM(ColorIndx(3),:),'box_on',1,...
%     'alphaval',5,'box_dodge',1,...
%     'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
%     'box_col_match',0,'box_dodge_amount',.5, 'dot_dodge_amount', 0.5,...
%     'line_width',3,'lwr_bnd',2);
% 
% raincloud_plot(pcaDYSK,'color',CM(ColorIndx(1),:),'box_on',1,...
%     'alphaval',5,'box_dodge',1,...
%     'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
%     'box_col_match',0,'box_dodge_amount',.75, 'dot_dodge_amount', 0.75,...
%     'line_width',3,'lwr_bnd',2);
% 
% pcaAMAN=Xpca(TWOpca,Y=='Amantadine');
% pcaCLZ=Xpca(TWOpca,Y=='Clozapine');
% pcaDYSK=Xpca(TWOpca,Y=='Dyskinesia');
% 
% figure
% raincloud_plot(pcaAMAN,'color',CM(ColorIndx(2),:),'box_on',1,...
%     'alphaval',5,'box_dodge',1,...
%     'box_dodge_amount',0.4, 'dot_dodge_amount', 0.4,...
%     'box_col_match',0,'box_dodge_amount',.25, 'dot_dodge_amount', 0.25,...
%     'line_width',3,'lwr_bnd',2);
% 
% raincloud_plot(pcaCLZ,'color',CM(ColorIndx(3),:),'box_on',1,...
%     'alphaval',5,'box_dodge',1,...
%     'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
%     'box_col_match',0,'box_dodge_amount',.5, 'dot_dodge_amount', 0.5,...
%     'line_width',3,'lwr_bnd',2);
% 
% raincloud_plot(pcaDYSK,'color',CM(ColorIndx(1),:),'box_on',1,...
%     'alphaval',5,'box_dodge',1,...
%     'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
%     'box_col_match',0,'box_dodge_amount',.75, 'dot_dodge_amount', 0.75,...
%     'line_width',3,'lwr_bnd',2);