% Function to get the Indexs of B according to Sorting in A sorted as C:
% such that C(1:end)==A(indxA)==B(indxB) as long as they are the same pair vectors
% if not-> indx's & C are empty
% Input
%   A:  [EXPIDs;Condition]  Size: Nx2 Categorical Vector
%   B:  [EXPIDs;Condition]  Size: Nx2 Categorical Vector
% Output
%   indx
function [C,indxA,indxB]=getindxofAinB(A,B)
% Setup
indxA=[];
indxB=[];
C=[];
if size(A,1)==size(B,1) && size(A,2)==size(B,2)
    disp('>>Vector Sizes: allright!')
    [C,indxA,indxB] = intersect(A,B,'rows');
    % Check Sorted Intersection:
    if size(A,1)==size(C,1) && size(A,2)==size(C,2)
        disp('>>Vector Sorting: allright!')
    else
        disp('>>Vector Sorting: ERROR: check missing or repeated experiments!')
        indxA=[];
        indxB=[];
        C=[];
    end
    
else
    disp('>>Vector Sizes: ERROR: check table mergeing')
end


