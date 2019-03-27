% Funtion To Deconvolve pair of Funtions
% x=d*r
% Input
%   x:      Signal
%   r:      Response Signal
%   lambda: Sparse Parameter
% Output
%   d:      Driver Signal
%   x_est:  Estimated Sparse Signal
%   MSEE:   Mean Squared Error
function [d,x_est,MSEE]=magic_sparse_deconvolution(x,r,lambda)
        
L=length(r);
Np=numel(x);
% To avoid probles with the size of the response
if L>Np
    r=r(1:Np);
    L=Np;    
end
% Build Response Matrix
Cols = [r';zeros(abs(length(x)-L),1)];
Rows = [r(1),zeros(1,length(x)-1)];
R = toeplitz(Cols,Rows);    
% Deconvolution
d = SpaRSAx(x',R,lambda);       % Driver 
% Estimation
x_est = R*d;                    % Sparse Estimation
MSEE=mean((x-x_est').^2);