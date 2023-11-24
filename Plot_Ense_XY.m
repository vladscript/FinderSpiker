XY=XY_selectedClean(Index_Ensemble,:);
C=fNet.Colors;
figure;
scatter(XY(:,1),XY(:,2),30*ones(size(XY(:,2))),C(Index_Ensemble,:)./255,'filled')
Ax=gca;
Ax.YDir="reverse";
Ax.XLim=[0,1.1*max(XY(:,1))];
Ax.YLim=[0,1.1*max(XY(:,2))];