% Model Multiclass SVM Classification binary complete
% Input:
%   Xpca:               Data matrix
%   ONEpca:             Column i-th of Xpca
%   TWOpca:             Column j-th of Xpca
%   Y:                  labels
%   template:           Kernel SVM template
%   labelconditions:    Cell of strings 
%   Nsim:               N models
% Output
%   PosterioirRegion:   Average probability regions
%   Ymany:              Nsim predictions
function [PosteriorRegion,Ymany]=svmmulticlass(Xpca,ONEpca,TWOpca,Y,template,...
    labelconditions,Nsim,x1Grid,x2Grid)
Ymany=[];
ProbSups=[];
for n=1:Nsim
    fprintf('>>Building model %i: ',n)
    Mdl = fitcecoc(Xpca(:,[ONEpca,TWOpca]),Y, 'Learners', template, 'Coding',...
            'binarycomplete', 'PredictorNames', {['PC_',num2str(ONEpca)],['PC_',num2str(TWOpca)]}, ...
            'ResponseName', 'Y','FitPosterior',1, ...
            'ClassNames', labelconditions);
    
    % maybe average model parameters ... (!)
        
    % Estimation of the labels:    
    [Yhat,~,~,~]=resubPredict(Mdl);
    % [Yhat,~,~,pPosterior]=resubPredict(Mdl);
    % possible use: [Xroc,Yroc,Troc,AUCroc] = perfcurve(Yhat,pPosterior,PosClass);
    Ymany=[Ymany,Yhat];
    fprintf('done.\n')
    fprintf('>>Predicting Regions: ... ')
    [~,~,~,PosteriorRegion] = predict(Mdl,[x1Grid(:),x2Grid(:)]);
    if n==1;
        ProbSups=PosteriorRegion;
    else
        ProbSups=ProbSups+PosteriorRegion;
    end
    fprintf('done.\n')
end
PosteriorRegion = ProbSups/Nsim;