% Function to plto conout of probabilities
function plot_multiclass(x1Grid,x2Grid,PosteriorRegion,CMprobs,Xpca,ONEpca,...
        TWOpca,Y,CM,ColorIndx,BestKernel,AccuracySVM,Yclass)
fprintf('>>Plotting regions: ... ')
MissClassifed=find(Yclass~=Y);
figure;
contourf(x1Grid,x2Grid,...
        reshape(max(PosteriorRegion,[],2),size(x1Grid,1),size(x1Grid,2)));
Contorno=gca;
colormap(Contorno,CMprobs);
h = colorbar;
h.YLabel.String = 'Maximum posterior';
h.YLabel.FontSize = 8;
hold on
% SCATTER PLOT  WITH DATA POINTS
gh = gscatter(Xpca(:,ONEpca),Xpca(:,TWOpca),Y,'krb','ooo',12);
gscatter(Xpca(MissClassifed,ONEpca),Xpca(MissClassifed,TWOpca),Y(MissClassifed),'kkk','xxx',12);
for n=1:numel(ColorIndx)
    gh(n).MarkerEdgeColor=[0,0,0];
    gh(n).MarkerFaceColor=CM(ColorIndx(n),:);
    gh(n).LineWidth=1;
end

title(['SCV kernel ',BestKernel,' LOOCV ACC=',num2str(round(1000*AccuracySVM)/10),'% ']);
xlabel(['PC ',num2str(ONEpca)]);
ylabel(['PC ',num2str(TWOpca)]);
% axis([xMin(1),xMax(1),xMin(2),xMax(2)])
axis([min(x1Grid(:)),max(x1Grid(:)),min(x2Grid(:)),max(x2Grid(:))]);
legend(gh,'Location','NorthWest');
hold off
fprintf(' done.\n')