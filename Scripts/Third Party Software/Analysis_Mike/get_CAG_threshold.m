% Input
%   raster: frames x cells (heritage from JP)
%   pval
% Output
%   Ci
% Example
function th=get_CAG_threshold(raster,alphaC)
% All in Same
Ncells=size(raster,1);
cluster_index=ones(Ncells,1);
th=testCoactivityGroup(raster,cluster_index,alphaC);