% Raster Ensembles Using Images FAST
% Input
%   R_Analysis: Raster Condition Structure
%           Raster: Frames x Cells
%   labels_frames
%   ColorState: Ensemble Colors
% Output
% Image with colored ensembles
function ImageEnsembles(R_Analysis,plothebb,varargin)
%% Setup
if ~isempty(varargin)
    ColorState=varargin{1}; % Alread Colored Ensembles
else
    TotalNG=R_Analysis.Clustering.TotalStates;
    NC=1;
    NGroups{1}=TotalNG;
    ColorState=colormapensembles(TotalNG,NC,NGroups);
end

R=R_Analysis.Data.Data; % Frams x Cells
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
%% Make Images
Rimage=R'; % Cells x Frames
CAG=sum(R');
CAGimage=zeros(max(CAG),size(R,1));
activityframes=find(CAG);
for f=1:numel(activityframes)
    CAGimage(CAG(activityframes(f)),activityframes(f))=1;
end
%% Colored Ensembles
for n=1:Nens
    % Find cells, frames and Nens+2;
    frames_ensemble=active_Frames( label_Cluster==n );
    Rimage(:,frames_ensemble)=R(frames_ensemble,:)'.*(1+n);
    CAGimage(:,frames_ensemble)=CAGimage(:,frames_ensemble).*(1+n);
end
%% Create Figure
EnsemFig=figure; % RASTER
EnsemFig.Name='Neural Ensembles';
axraster=subplot(3,1,[1,2]);
imagesc(Rimage);
axcag=subplot(3,1,3);
% plot(CAG,'k'); hold on;
imagesc(CAGimage); % hold off;
colormap(axraster,ColorEnsembles)
colormap(axcag,ColorEnsembles)
axraster.YDir='normal';
axcag.YDir='normal';
% figure; % Hebbian Ensembles
Ensembles_Transitions(1/60,label_Cluster,active_Frames,ColorEnsembles(3:end,:),plothebb); % ---> save
linkaxes([axraster,axcag],'x');
