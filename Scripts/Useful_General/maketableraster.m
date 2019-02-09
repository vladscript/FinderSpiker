% function that builds a table to save ti as CSV
% from several data cells of Raster Features
% Input
%   RASTER_NAME:        cell with ExpIDs
%   RASTER_FEATURES:    cell  of rastrfeatures 
%   Names_Conditions:   cell of Condition Names
function TableFeatures=maketableraster(RASTER_NAMES,RASTER_FEATURES,Names_Conditions)
NC=numel(Names_Conditions);
TableFeatures=table;
for c=1:NC
    Nexps=numel(RASTER_NAMES{c});
    ActualName={};
    for n=1:Nexps
        ActualName{n,1}=Names_Conditions{c};
    end 
    TableFeatures=[TableFeatures;table(ActualName,RASTER_NAMES{c},RASTER_FEATURES{c})];
end
