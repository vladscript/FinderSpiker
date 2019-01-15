% Raster Ensembles Using Images FAST
% Input
%   R_Analysis: Raster Condition Structure
%           RasterFrames x Cells
%   labels_frames
%   ColorState: Ensemble Colors
% Output
% Image with colored ensembles
function plot_ensambles_fast(R_Analysis,ColorState)
%% Setup
R=R_Analysis.Data.Data;
label_Cluster=R_Analysis.Clustering.VectorStateIndex;
CAGth=R_Analysis.Peaks.Threshold;
active_Frames=find(sum(R,2)>=CAGth);
[NensPlus,~]=size(ColorState);
Nens=NensPlus-1;
% ScaleValues=0:Nens+1; % background and non ensemble activity
ColorEnsembles=[1,1,1;0,0,0;ColorState(1:Nens,:)];
%               white; black; EnsembleColors]
% Initialize Oouput
% Rimage=zeros(size(R));
% [TotalFrames,~]=size(R);
Rimage=R'; % Cells x Frames
CAGimage=zeros(max(sum(Rimage,2)))
for n=1:Nens
    % Find cells, frames and Nens+2;
    frames_ensemble=active_Frames( label_Cluster==n );
    Rimage(:,frames_ensemble)=R(frames_ensemble,:)'.*(1+n);
end
%% Create Figure
figure; 
axraster=subplot(3,1,[1,2]);
imagesc(Rimage);
colormap(axraster,ColorEnsembles)
axraster.YDir='normal';
axCAG=subplot(3,1,3);

