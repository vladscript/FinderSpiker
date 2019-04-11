% Script that gets the Network from Whole Raster
% Input:
%   R_Condition: Cell of Rasters
%   MetaDataColocaliation: List of [+] & [-] Cells
% Outputs:
% Colored Raster
% CS files to get the Gephi Network
% Update .mat File
C=size(R_Condition,2);
if exist('MetaDataColocaliation')
    aremerged=true;   % Are there already colocated cells
    PositiveCells=MetaDataColocaliation.PositiveCells;
    NegativeCells=MetaDataColocaliation.NegativeCells;
else
    aremerged=false;  % Are there already colocated cells
end

for c=1:C
    Raster=R_Condition{c};
    Ncells=size(Raster,1);
    % TOTAL NETWORK
    AdjacencyMatrix=GetAdjacencyMatrix(Raster);
    MaxSynLinks=max(max(AdjacencyMatrix));
    AdjacencyMatrix=AdjacencyMatrix./MaxSynLinks; % NORMALIZED
    % Get Source, Target and Weigth for eache link NO SORTED
    SOURCE=[];
    TARGET=[];
    WEIGHT=[];
    for i=1:Ncells-1
        for j=i+1:Ncells
            if AdjacencyMatrix(i,j)>0
                SOURCE=[SOURCE;i];
                TARGET=[TARGET;j];
                WEIGHT=[WEIGHT;AdjacencyMatrix(i,j)];
            end
        end
    end
    if isempty(SOURCE)
        SOURCE=0;
        TARGET=0;
        WEIGHT=0;
        disp('>>NO NETWORK FOUND')
    end
    
    if aremerged
    else
    end
end
    