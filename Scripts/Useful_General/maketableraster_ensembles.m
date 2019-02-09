% function that builds a table to save ti as CSV
% from several data cells of Raster Features
% Input
%   RASTER_NAME:        cell with ExpIDs
%   RASTER_FEATURES:    cell  of rastrfeatures 
%   Names_Conditions:   cell of Condition Names
function TableFeatures=maketableraster_ensembles(RASTER_NAMES,RASTER_FEATURES,NensemblesCondition,Names_Conditions)
NC=numel(Names_Conditions);
TableFeatures=table;
for c=1:NC
    TotalNames=[];
    Nexps=numel(RASTER_NAMES{c});
    Nfeats=size(RASTER_FEATURES{c},1);
    for e=1:Nexps
        Nens=NensemblesCondition{c}(e);
        ActualName={};
        for k=1:Nens
            ActualName{k,1}=RASTER_NAMES{c}(e);
        end
        TotalNames=[TotalNames;ActualName];
    end
    ActualCondition={};
    for k=1:Nfeats
        ActualCondition{k,1}=Names_Conditions{c};
    end
    
    TableFeatures=[TableFeatures;table(ActualCondition,TotalNames,RASTER_FEATURES{c})];
end
