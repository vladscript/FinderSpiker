% Get Adjacency with the connected nodes

function [A IdxConnected IdxNotConnected] = OnlyConnected(Adjacency)
    Links=sum(Adjacency);
    IdxConnected=find(Links>0);
    A=Adjacency(IdxConnected,IdxConnected);
    A=double(A>0);
    IdxNotConnected=find(Links==0);
end