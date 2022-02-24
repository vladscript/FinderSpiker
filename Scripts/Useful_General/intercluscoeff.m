%% Interclsuter Coefficient (ICC) & Effective N (Neff)
% As defined in AArts et al. 2021
% 
% Input
%   X: numerical input  data
%   Y: category of input  data
% Ouput
%  ICC: VarB/[VarB + VarW], where:
%       VarB: variance between research objects
%       VarW: the sum of the variance within research objects
%
% ICC=0-> all observations obtained from a research object are independent
% ICC=1-> all observations obtained from the same research objects are equal 
% and therefore convey the very same information.
% 
% Example:
% 9 samples of 3 measures of 3 samples (A. B and C
% X=[1,2,3,4,5,6,7,8,9];
% Y={'A','A','A','B','B','B','C','C','C'};
% ICC==intercluscoeff(X,Y);
% 
% Alternatively Y can be:
%   Y={'A','A','A','B','B','B','C','C','C'};
%   Y=categorical({'A','A','A','B','B','B','C','C','C'})
%   Y=['A','A','A','B','B','B','C','C','C'];
% 
function [ICC,Neff]=intercluscoeff(X,Y)
% Settings ###############################################
if ischar(Y)
    Ystr=Y; Y=cell(numel(Y),1);
    for i=1:numel(Y)
        Y{i}=Ystr(i);
    end
end
Yunique=unique(Y,'stable');
N=numel(X);
Nclusters=numel(Yunique);
meanClusts=zeros(Nclusters,1);

tabulate(Y);
% Main Loop ###############################################
for n=1:Nclusters
    if iscategorical(Y)
        meanClusts(n)=mean(X(Y==Yunique(n)));
        Nsamp(n)=numel(X(Y==Yunique(n)));
        varClus(n)=var(X(Y==Yunique(n)));
    elseif iscell(Y)
        meanClusts(n)= mean(X(ismember(Y,Yunique{n}))); 
        Nsamp(n)=numel(X(ismember(Y,Yunique{n})));
        varClus(n)= var(X(ismember(Y,Yunique{n}))); 
    else
        disp('>Check data in Y')
        meanClusts=NaN;
        Nsamp=0;varClus=0;
    end
    
end
VarB=var(meanClusts);
VarW=sum(varClus);
ICC=VarB/(VarB+VarW);
Neff=(1./(1+(Nsamp)*ICC))*N;
% #########################################################