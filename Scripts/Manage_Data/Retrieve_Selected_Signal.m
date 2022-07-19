% Function to Retrieve Signals From Selected Frames of the Experiment
% Input
%   Onsets:         Cell of Sample Onsets per Condition
%   R_Condition:    Selected Raster per Condition (cell)
%   X:              Original Signals (cell)
%   XY:             Original Set of Coordinates (Matrix)
%   XY_selected:    Selected Coordinates (Matrix)
% Such that:         XY_selected = XY(IndexSorted,:)
% Output
%   Indexes:    SELECTED Original Indexes
%   Xsel:       Concatenated Selected Signals
function [X_sel,IndexSorted]=Retrieve_Selected_Signal(Onsets,R_Condition,X,XY_selected,XY)
% Indexes of Original that were Selected
[~,~,IndexSorted]=intersect(XY_selected,XY,'rows','stable');
% Loook for repeated ones:
IndxRepeated=repxychekcker(XY_selected);

[NV,NC]=size(X);
X_sel=cell(1,NC);
for n=1:NC
    Start=Onsets{n};
    [~,End]=size(R_Condition{n});
    Xmat=[];
    for v=1:NV
        if ~isempty(X{v,n})
            Xmat=[Xmat,X{v,n}];
        end
    end
    X_sel{n}=Xmat(IndexSorted,Start:Start+End-1);
end
end