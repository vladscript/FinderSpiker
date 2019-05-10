%% script gets mean PDF of Distance:
% between or among different conditions of
% rasters:
% A vector is fromed where Active Cells makes 1
% Input
%   R_Condition: Cell of rasters
%   Distance (from matlab): 'hamming' (default), 'euclidean', 'cosine'
% Output
% MatDist: Matrix Of Dsitacne among Conditions:
%   MatDis(1,2) is the distance bewtween COndition 1 and 2
function MatDist=distance_rasters(R_Condition,varargin)
if ~isempty(varargin)
    NameDistance=varargin{1};
else
    NameDistance='hamming';
end
fprintf('Distance to measure: %s\n',NameDistance);
% Setup
NC=size(R_Condition,2);
Ncells=size(R_Condition{1},1);
% Get Binary vector 
if NC>1
    AC=zeros(Ncells,NC);
    for c=1:NC
        R=R_Condition{c};
        AC(sum(R,2)>0,c)=1; 
    end
    % calculate Distance Matrix
    MatDist=squareform(pdist(AC',NameDistance));
else
    disp('>>Not Useful for only 1 Condition');
    MatDist=0;
end