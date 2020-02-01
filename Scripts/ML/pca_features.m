%% Get Data

% X=copy data from ALLDATA.xlsx

Y(Y=='DyskinesiaA')='Dyskinesia';
Y(Y=='DyskinesiaC')='Dyskinesia';

Y = removecats(Y); % removes empty categories

%% PCA
% function Xpca=pca_features(X,varLevel)
[coefs,score,latent] = pca(X,'algorithm','als');
% Rows of score correspond to observations, 
% and columns to components.
VarExplained=cumsum(latent)./sum(latent);
Npcs=find(VarExplained>=0.95,1);
Xpca=coefs*score'; % Features c observations
[Nfeatures,Nobser]=size(Xpca);

%% BOXPLOTS

[CM,ColorIndx]=Color_Selector({'Dyskinesia','+Amantadine','+Clozapine'});

figure
for i=1:16
    pcaAMAN=Xpca(i,Y=='Amantadine');
    pcaCLZ=Xpca(i,Y=='Clozapine');
    pcaDYSK=Xpca(i,Y=='Dyskinesia');
    ksdensity(pcaAMAN,linspace(min(pcaAMAN),max(pcaAMAN),100)); hold on;
    ksdensity(pcaCLZ,linspace(min(pcaCLZ),max(pcaCLZ),100));
    ksdensity(pcaDYSK,linspace(min(pcaDYSK),max(pcaDYSK),100)); hold off;
    axis tight; grid on;
    legend('Amantadine','Clozapine','Dyskinesia');
    %     columns2boxplot(pcaDYSK',pcaAMAN',pcaCLZ',{'dysk','aman','clz'})
    figure
    raincloud_plot(pcaAMAN,'color',CM(ColorIndx(2),:),'box_on',1,...
        'alphaval',5,'box_dodge',1,...
        'box_dodge_amount',0.4, 'dot_dodge_amount', 0.4,...
        'box_col_match',0,'box_dodge_amount',.25, 'dot_dodge_amount', 0.25,...
        'line_width',3,'lwr_bnd',2);

    raincloud_plot(pcaCLZ,'color',CM(ColorIndx(3),:),'box_on',1,...
        'alphaval',5,'box_dodge',1,...
        'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
        'box_col_match',0,'box_dodge_amount',.5, 'dot_dodge_amount', 0.5,...
        'line_width',3,'lwr_bnd',2);

    raincloud_plot(pcaDYSK,'color',CM(ColorIndx(1),:),'box_on',1,...
        'alphaval',5,'box_dodge',1,...
        'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
        'box_col_match',0,'box_dodge_amount',.75, 'dot_dodge_amount', 0.75,...
        'line_width',3,'lwr_bnd',2);
    pause
end



%% Classification 

template = templateSVM('KernelFunction', 'linear', 'PolynomialOrder', [],...
    'KernelScale', 'auto', 'BoxConstraint', 1, 'Standardize', 1,'Verbose',1);

aux=0; PCAandACC=[];
for i=1:Nfeatures-1
    for j=i+1:Nfeatures
        % SVM
        Mdl = fitcecoc(Xpca([i,j],:)',Y, 'Learners', template, 'Coding',...
        'onevsone', 'PredictorNames', {'Var1_1' 'Var1_9'}, ...
        'ResponseName', 'Y','FitPosterior',1, ...
        'ClassNames', categorical({'Amantadine' 'Clozapine' 'Dyskinesia'}));
        %[Yhat,~,~,pPosterior]=resubPredict(Mdl);
        % Evlauate 
        partitionedModel = crossval(Mdl, 'KFold', 5);
        validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
        % Save results
        PCAandACC=[PCAandACC;i,j,validationAccuracy];
        fprintf(repmat('*',20,1));
        fprintf('\n>>Model %i, PCA components: %i vs %i\n',aux,i,j);
        fprintf(repmat('*',20,1));
        aux=aux+1;
    end
end
% GET BEST PCAs
[~,nModel]=max(PCAandACC(:,3));
ONEpca=PCAandACC(nModel(1),1);
TWOpca=PCAandACC(nModel(1),2);

Mdl = fitcecoc(Xpca([ONEpca,TWOpca],:)',Y, 'Learners', template, 'Coding',...
        'onevsone', 'PredictorNames', {'Var1_1' 'Var1_9'}, ...
        'ResponseName', 'Y','FitPosterior',1, ...
        'ClassNames', categorical({'Amantadine' 'Clozapine' 'Dyskinesia'}));
    
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
plotroc(Targetss',Outputss')
plotconfusion(Targetss',Outputss')

%% RAINCLOUDS of THE BEST PCAs
pcaAMAN=Xpca(ONEpca,Y=='Amantadine');
pcaCLZ=Xpca(ONEpca,Y=='Clozapine');
pcaDYSK=Xpca(ONEpca,Y=='Dyskinesia');

figure
raincloud_plot(pcaAMAN,'color',CM(ColorIndx(2),:),'box_on',1,...
    'alphaval',5,'box_dodge',1,...
    'box_dodge_amount',0.4, 'dot_dodge_amount', 0.4,...
    'box_col_match',0,'box_dodge_amount',.25, 'dot_dodge_amount', 0.25,...
    'line_width',3,'lwr_bnd',2);

raincloud_plot(pcaCLZ,'color',CM(ColorIndx(3),:),'box_on',1,...
    'alphaval',5,'box_dodge',1,...
    'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
    'box_col_match',0,'box_dodge_amount',.5, 'dot_dodge_amount', 0.5,...
    'line_width',3,'lwr_bnd',2);

raincloud_plot(pcaDYSK,'color',CM(ColorIndx(1),:),'box_on',1,...
    'alphaval',5,'box_dodge',1,...
    'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
    'box_col_match',0,'box_dodge_amount',.75, 'dot_dodge_amount', 0.75,...
    'line_width',3,'lwr_bnd',2);

pcaAMAN=Xpca(TWOpca,Y=='Amantadine');
pcaCLZ=Xpca(TWOpca,Y=='Clozapine');
pcaDYSK=Xpca(TWOpca,Y=='Dyskinesia');

figure
raincloud_plot(pcaAMAN,'color',CM(ColorIndx(2),:),'box_on',1,...
    'alphaval',5,'box_dodge',1,...
    'box_dodge_amount',0.4, 'dot_dodge_amount', 0.4,...
    'box_col_match',0,'box_dodge_amount',.25, 'dot_dodge_amount', 0.25,...
    'line_width',3,'lwr_bnd',2);

raincloud_plot(pcaCLZ,'color',CM(ColorIndx(3),:),'box_on',1,...
    'alphaval',5,'box_dodge',1,...
    'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
    'box_col_match',0,'box_dodge_amount',.5, 'dot_dodge_amount', 0.5,...
    'line_width',3,'lwr_bnd',2);

raincloud_plot(pcaDYSK,'color',CM(ColorIndx(1),:),'box_on',1,...
    'alphaval',5,'box_dodge',1,...
    'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
    'box_col_match',0,'box_dodge_amount',.75, 'dot_dodge_amount', 0.75,...
    'line_width',3,'lwr_bnd',2);

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

title 'Best 2-classifying PCA components and Maximum Posterior';
xlabel([num2str(ONEpca),'th PCA component']);
ylabel([num2str(TWOpca),'th PCA component']);
axis tight
legend(gh,'Location','NorthWest')
hold off
%% GSCATTERS per CONDITION one-versus-all (OVA)
figure
gscatter(Xpca(ONEpca,:),Xpca(TWOpca,:),Y,CM([2,3,1],:),'dso',10);
title('{\bf Best calssifyng PCA components}');
xlabel([num2str(ONEpca),' th PCA component']);
ylabel([num2str(TWOpca),' th PCA component']);
legend('Location','Northwest'); 
axis tight

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





%% Other STUFF
partitionedModel = crossval(trainedClassifier, 'KFold', 5);
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

[label,~,~,Posterior] = resubPredict(Mdl,'Verbose',1);
Mdl.BinaryLoss

idx = randsample(size(Xpca(1,:)',1),10,1);
Mdl.ClassNames
table(Y(idx),label(idx),Posterior(idx,:),...
    'VariableNames',{'TrueLabel','PredLabel','Posterior'})

xMax = max(Xpca(1,:)');
xMin = min(Xpca(1,:)');

x1Pts = linspace(xMin(1),xMax(1));
x2Pts = linspace(xMin(1),xMax(1));
[x1Grid,x2Grid] = meshgrid(x1Pts,x2Pts);

[~,~,~,PosteriorRegion] = predict(Mdl,[x1Grid(:);x2Grid(:)]);



ifeat=1;
% FeatureName=unique(Y,'legacy');

figure;
contourf(x1Grid,x2Grid,...
        reshape(max(PosteriorRegion,[],2),size(x1Grid,1),size(x1Grid,2)));
h = colorbar;
h.YLabel.String = 'Maximum posterior';
h.YLabel.FontSize = 15;
hold on
gh = gscatter(Xpca(1,:)',Xpca(1,:)',Y,'grkgb','do.ox',10);
gh(1).Color=CM(2,:);
gh(1).LineWidth=2;
gh(2).Color=CM(3,:);
gh(2).LineWidth=2;
gh(5).Color=CM(1,:);
gh(5).LineWidth=2;
% gh(2).LineWidth = 2;
% gh(3).LineWidth = 2;

title('Neuronal Activity and Maximum Posterior');
xlabel '1st PCA component';
ylabel '1st PCA component';
axis tight
legend(gh([1,4,5]),'Location','NorthWest')
hold off


% %% Loop-Feature for Naive Bayes Classifier
% 
% for Npcs=1:3
%     Mdl=fitcnb(Xpca(1:Npcs,:)',Y,'DistributionNames','kernel');
%     [Yhat,~]=resubPredict(Mdl);
%     ErrorSelFeat(Npcs)=1-numel(find(Y==Yhat))/numel(Y);
%     [C,order]=confusionmat(Y,Yhat);
%     
%     disp(Npcs)
% end