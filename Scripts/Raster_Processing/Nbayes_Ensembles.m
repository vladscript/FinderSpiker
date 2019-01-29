% Function to evaluate Hierarchichal Clustering bye means of
% Naive Bayes Classifier
% Input
%   Rclust: Raster of active frames
%   frames_labels: Labels at each frame of the raster
% Output
%   Mdl: Classiffier
%   ECV: Error Cross Validation
function [Mdl,ECV]=Nbayes_Ensembles(Rclust,FRAMES_LABELS)
%% Setup
Nensembles=numel(unique(FRAMES_LABELS));
NeuralEnsemble=unique(FRAMES_LABELS);
% Labels
Y=categorical(FRAMES_LABELS);
% Features
X=Rclust';
% FeaturesOK=cell(Nensembles,1);
% for n=1:Nensembles
%     ActualEnsemble=NeuralEnsemble(n);
%     FeaturesOK{n}=find(sum(X(FRAMES_LABELS==ActualEnsemble,:)));
%     % LabelsNames{n}=['Ensemble_',num2str(ActualEnsemble)];
% end
%% Magic Classification
% Use fo kernel due to possible var(NeuralActivity)=0
disp('Ensembles Data:')
tabulate(Y);
disp('>>Training Classiffier:')
Mdl=fitcnb(X,Y,'DistributionNames','kernel');
disp('>>Classiffier Trained.')
disp('>>Validating....')
% 10-fold Cross Validation
[label,~]=resubPredict(Mdl);
ECV=1-numel(find(label==Y))/numel(Y);
fprintf('Cross-Validated Classification Error: %3.1f %%\n',100*ECV);
disp('>>Done.')

% Equivalent to do:
% CVMdl=crossval(Mdl);
% CVMdl.kfoldLoss;

% USEFUL:
% C=confusionmat(Y,label);

