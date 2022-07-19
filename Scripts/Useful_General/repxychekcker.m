% Functions that check repeated pair of coordinates
% Ouput
%   repetidas: indexes of repeated coordinates
% Input
%   Matrix of Coordinates [X,Y]
function IndxRepeated=repxychekcker(XY_selected)
IndxRepeated=[];
[~,IndexSortedUN]=unique(XY_selected,'rows','stable');
if numel(IndexSortedUN)<size(XY_selected,1)
    fprintf('\n>DUPLICATED COORDINATES\n')
    RepIndx=setdiff(1:size(XY_selected,1),IndexSortedUN);
    fprintf('Signal - X , Y \n');
    for n=1:numel(RepIndx)
    %         find(XY_selected==XY_selected(RepIndx(n),:)
        repetidas=find(sum(ismember(XY_selected,XY_selected(RepIndx(n),:)),2)==2);
        for i=1:size(repetidas,1)
            fprintf(' %i ,%3.1f, %3.1f\n',repetidas(i),XY_selected(repetidas(i),1),XY_selected(repetidas(i),2));
        end
        IndxRepeated=[IndxRepeated;repetidas];
    end
else
    fprintf('\n>ALL UNIQUE COORDINATES\n')
end