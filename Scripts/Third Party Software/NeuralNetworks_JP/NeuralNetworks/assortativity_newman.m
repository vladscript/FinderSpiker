%% assortativity_newman Assortativity Coefficient according to Newman
%This function calcule the Assortativity Coefficient that may explain how
%the relations are made in the network if it is between nodes of higher
%degrees or lower, or if the higher degree nodes are connected with the
%lower.


%lm is the adjancecy matrix

%r is the newman coefficient 

%If r > 0, the network is claimed to be assortative mixing; 
%while if r < 0, the network is called disassortative mixing.

function [r]=assortativity_newman(lm)
L=sum(lm)-transpose(diag(lm)); %contains a degree vector for each node
M=0.5*sum(L);
N=length(lm);
sum1=0;
sum2=0;
sum3=0;
for k=1:N-1;
    for j=k+1:N 
        if(lm(k,j)~=0)
     
            sum1=L(k)*L(j)+sum1;
            sum2=L(k)+L(j)+sum2;
            sum3=L(k)^2+L(j)^2+sum3;
            
        end
    end
end

a=(sum1*M^-1)-(0.5*M^-1*sum2)^2;
b=(.5*sum3*M^-1)-(0.5*M^-1*sum2)^2;

r=a/b;
end