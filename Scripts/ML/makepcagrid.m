% Make grid of values for SVM prediction
% 
% Input:
%   Xpca    matrix pf Principal Components
%   PC_i    i-th PC
%   PC_j    j-th PC
%   delta   Resolution
%   Npts (optional): default=100
% Output
%   
% 
function [x1Grid,x2Grid]=makepcagrid(Xpca,ONEpca,TWOpca,deltaPlus,varargin)
% N pts:
if isempty(varargin)
    Np=100;
else
    Np=varargin{1};
end
MaxFactorPlus=ones(1,2);
MinFactorPlus=MaxFactorPlus;
maxSigns=sign(max(Xpca(:,[ONEpca,TWOpca])));
minSigns=sign(min(Xpca(:,[ONEpca,TWOpca])));
MaxFactorPlus(maxSigns>0)=1+deltaPlus;
MaxFactorPlus(maxSigns<0)=1-deltaPlus;
MinFactorPlus(minSigns>0)=1-deltaPlus;
MinFactorPlus(minSigns<0)=1+deltaPlus;
xMax = max(Xpca(:,[ONEpca,TWOpca])).*MaxFactorPlus;
xMin = min(Xpca(:,[ONEpca,TWOpca])).*MinFactorPlus;
x1Pts = linspace(xMin(1),xMax(1),Np);
x2Pts = linspace(xMin(2),xMax(2),Np);
[x1Grid,x2Grid] = meshgrid(x1Pts,x2Pts);