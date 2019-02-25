%% assortativity_newman Assortativity Coefficient according to Newman
%This function calcule the Assortativity Coefficient that may explain how
%the relations are made in the network if it is between nodes of higher
%degrees or lower, or if the higher degree nodes are connected with the
%lower.


%A is the adjancecy matrix

%r is the newman coefficient 

%If r > 0, the network is claimed to be assortative mixing; 
%while if r < 0, the network is called disassortative mixing.


%Modified by Jesús E. Pérez-Ortega Nov-2012
function r=Assortativity(A)
L=sum(A)-transpose(diag(A)); %contains a degree vector for each node
M=0.5*sum(L);
N=length(A);
sum1=0;
sum2=0;
sum3=0;

B=triu(A,1);

for j=1:N
    idx=find(B(j,:));
    N2=length(idx);
    for k=1:N2;
        sum1=L(j)*L(idx(k))+sum1;
        sum2=L(j)+L(idx(k))+sum2;
        sum3=L(j)^2+L(idx(k))^2+sum3;
    end
end
a=(sum1*M^-1)-(0.5*M^-1*sum2)^2;
b=(.5*sum3*M^-1)-(0.5*M^-1*sum2)^2;

r=a/b;
end