% Scale free network
%
% This function generates binary undirected scale free network.
%
% SF = SF_JP(N, K)
%
% Inputs:
% N = number of nodes
% K = number of links (~ K links)
%
% Outputs:
% SF = binary undirected adjacent matrix of scale free network
%
% ..:: by Jesús E. Pérez-Ortega ::.. March-2013

function SF = ScaleFree_JP(N,K)
    SF=eye(N);
    m=round(min(roots([1 -N K]))); % What if complex roots?
    P=ones(m,1)/m;
    for i=1:N-m
        node=m+i;
        for j=1:m
            idxP=find(~SF(node,1:node-1));
            Pnodes=P(idxP)/sum(P(idxP));
            idx_n2=1+sum(rand>cumsum(Pnodes));
            node2=idxP(idx_n2);
            SF(node,node2)=1;
            SF(node2,node)=1;
            P=sum(SF)/sum(sum(SF));
        end
    end
    SF=SF-eye(N);
end