% Generate Ring Regular Network
%
% This function generates binary undirected network with ring regular
% connections.
%
% Reg = RegNet(N,K)
%
% Inputs:
% N = number of nodes
% K = number of links
%
% Outputs:
% Reg = binary undirected adjacent matrix of ring regular network
%
% ..:: by Jesús E. Pérez-Ortega ::.. March-2013

function Reg = RegNet_JP(N,K)
    k=2*K/N;
    Reg=zeros(N);
    km=round(ceil(k)/2);
    for i=1:N
        n_b=i-[km:-1:1];
        n_a=i+[1:km];
        
        initial=find(n_b<1);
        if initial
            n_init=length(initial);
            n_b(1:n_init)=[N:-1:(N-n_init+1)];
        end
        
        final=find(n_a>N);
        if final
            n_final=length(final);
            n_a(end:-1:(end-n_final+1))=[1:n_final];
        end
        
        Reg(i,[n_b, n_a])=1;
        Reg([n_b, n_a],i)=1;
    end
    
    % Remove excess connections
    n_remove=(sum(sum(Reg))-ceil(k*N))/2;
    [i j]=find(triu(Reg));
    idx=round(rand(1,n_remove)*length(i));
    for l=1:n_remove
        Reg(i(idx(l)),j(idx(l)))=0;
        Reg(j(idx(l)),i(idx(l)))=0;
    end
end