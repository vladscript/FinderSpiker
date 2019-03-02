% function Xpca=pca_features(X,varLevel)
[coefs,score,latent] = pca(X,'algorithm','als');
% Rows of score correspond to observations, 
% and columns to components.


VarExplained=cumsum(latent)./sum(latent);
Npcs=find(VarExplained>=0.95,1);
Xpca=coefs*score';

for Npcs=1:3
    Mdl=fitcnb(Xpca(1:Npcs,:)',Y,'DistributionNames','kernel');
    [Yhat,~]=resubPredict(Mdl);
    ErrorSelFeat(Npcs)=1-numel(find(Y==Yhat))/numel(Y);
    [C,order]=confusionmat(Y,Yhat)
    disp(Npcs)
end