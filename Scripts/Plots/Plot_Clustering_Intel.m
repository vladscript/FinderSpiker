%% TEST ILLLSUTRATE
Nensambles=6;
frame_ensembles=cluster(Tree,'maxclust',Nensambles); % START
% Checking Redundant Ensembles
NeuroVectors=zeros(20,Nensambles);
for nn=1:Nensambles
    EnsembledNeurons{nn}=find(sum(Rclust(:,frame_ensembles==nn),2));
    NeuroVectors(EnsembledNeurons{nn},nn)=1;
end
DistaVec=squareform(pdist(NeuroVectors',SimMethod));
Plot_Raster_Ensembles(NeuroVectors);
drawnow;
Rfig=gcf;
Raxis=Rfig.Children(end);
Raxis.YTick=1:20;
Raxis.XLim=[0.5,5.5];
SimVectors=1-DistaVec

%% SCRIPT
% Plot Clustering Intel
% Get CAG threshold
%   Rclust
% Plot Simm Matrix
% Get Linkage Method
% Plot Simm Matrix
% Sort Raster by Ensembles
% Load_Default_Clustering;
Distance=squareform(pdist(Rclust',SimMethod));
Sim=1-Distance;
figure; imagesc(Sim); % SIMM MATRIX
LinkageMethod=HBestTree_JPplus(Sim);    % Output
Tree=linkage(squareform(Distance,'tovector'),LinkageMethod);
figure;
[H,T,outperm]=dendrogram(Tree,0,'Orientation','right',...
    'ColorThreshold',0.25);

figure; 
imagesc(Sim(T(outperm(end:-1:1)),T(outperm(end:-1:1)))); % SIMM MATRIX SORTED
