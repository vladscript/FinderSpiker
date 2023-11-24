% RUN AFTER >>Ensemble_Sorting

XY=XY_selectedClean(Index_Ensemble,:);
C=fNet.Colors;
ColorsInd=C(Index_Ensemble,:)./255;
for i=1:size(ColorsInd,1)
    actcl=ColorsInd(i,:);
    axx=find(ismember(ColorState,actcl,'rows'));
    if ~isempty(axx)
        Iensamble(i)=axx;
    else
        Iensamble(i)=0;
    end
end


figure;
scatter(XY(:,1),XY(:,2),30*ones(size(XY(:,2))),ColorsInd,'filled')
Ax=gca;
Ax.YDir="reverse";
Ax.XLim=[0,1.1*max(XY(:,1))];
Ax.YLim=[0,1.1*max(XY(:,2))];

% Generara Tabla CSV con coordenadas y grupo
Txyens=table(XY(:,1),XY(:,2),Iensamble');
Txyens.Properties.VariableNames={'X','Y','EnsembleS'};
SaveDirectory=uigetdir(pwd,'Select Folder Destinty for coordinates table:');
if SaveDirectory==0
    fprintf('\n>Table unsaved\n')
else
    writetable(Txyens,[SaveDirectory,filesep,Experiment,'_XY_Ensembles.csv']);
    fprintf('\n>Table saved\n')
    fprintf('<a href="matlab:dos(''explorer.exe /e, %s, &'')">See coordinates here</a>\n',SaveDirectory);
end
