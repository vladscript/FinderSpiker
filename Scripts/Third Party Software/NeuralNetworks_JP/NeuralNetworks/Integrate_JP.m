% Integrate
% Compute the integrate of data.
%
% [I] = Integrate_JP(X)
%
% Inputs
% X = data as vector of n elements.
%
% Output
% I = integral of the data.
%
% ..:: by Jesús E. Pérez-Ortega ::.. May-2012

function [I] = Integrate_JP(X)

n=length(X);
I=zeros(1,n);

for i=1:n
    I(i)=sum(X(1:i));
end