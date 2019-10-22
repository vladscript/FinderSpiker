% Sparse Convolution: x_sparse=R*d
% Input
%   d: drive vector
%   r: response vector
% Output
%   x_sparse: drive vector
function x_sparse=sparse_convolution(d,r)
d=makerowvector(d)';
r=makerowvector(r)';
L=numel(r);
Np=numel(d);
% To avoid problems with the size of the response
if L>Np
    r=r(1:Np);
    L=numel(r);    
end
% Build Response Matrix
Cols = [r;zeros(abs(Np-L),1)];
Rows = [r(1),zeros(1,Np-1)];
R = toeplitz(Cols,Rows); 
x_sparse = R*d;
end