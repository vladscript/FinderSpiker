%% Get Connection
function [Connected Disconnected NumCon NumDis]=Get_Connection(Adjacency)
    Disconnected=find(sum(Adjacency,1)==0);
    NumDis=length(Disconnected);
    Connected=find(sum(Adjacency,1)>0);
    NumCon=length(Connected);
end