function ShowmeyourCoordinates(XY,r,isSIGNAL)
%% Show Coordinates
% ALL
figure;
scatter(XY(:,1),XY(:,2),r.*r*pi,'LineWidth',1.5)
Ncells=numel(r);
Indexes=1:Ncells;
hold on
scatter(XY(isSIGNAL,1),XY(isSIGNAL,2),r(isSIGNAL).*r(isSIGNAL)*pi,...
    repmat([0.2,0.3,0.5],numel(isSIGNAL),1),'filled','LineWidth',1.5)
for n=1:Ncells
    text(XY(n,1)+round(r(n)/2),XY(n,2),num2str(n))
end
grid on
ActualFig=gcf;
ActualFig.Children.YDir='reverse';
% SELECTED