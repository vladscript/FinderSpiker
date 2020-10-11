%% Plot Venn Diagram
%  Input
%   ActivityGroups
% Output
%   Plot
function plotvennCellActive(ActiveCells,Names_Conditions)
% 3 main SEPARATED groups:
KeepDoPlot=true;
if max(union(ActiveCells.AC_1,ActiveCells.AC_2))==max(ActiveCells.CellsIndx)
    fprintf('>>OK')
else
    fprintf('Something Wrong.\n');
    KeepDoPlot=false;
end


if KeepDoPlot
    %% Areas & Intersections
    TotalArea=300;
    Acond1=numel(ActiveCells.AC_1)/numel(ActiveCells.CellsIndx)*TotalArea;
    Acond2=numel(ActiveCells.AC_2)/numel(ActiveCells.CellsIndx)*TotalArea;
    % Select Color 
    [CM,ColorIndx]=Color_Selector(Names_Conditions);
    Inter=numel(intersect(ActiveCells.AC_1,ActiveCells.AC_2))/numel(ActiveCells.CellsIndx)*TotalArea;
    figure%, axis equal, axis off
    venn([Acond1,Acond2],Inter,'FaceColor',{CM(ColorIndx(1),:),CM(ColorIndx(2),:)},...
        'FaceAlpha',{1,0.6},'EdgeColor','black')
    legend(Names_Conditions);
    axis image
    
end