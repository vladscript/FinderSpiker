% Get Adjacency
function [Adjacency ID W_Th WAdjacency]= Get_Adjacency(Data,P_Th,W_Th,W_Th_Auto)
    % Inputs
    Data=double(Data);
    C=size(Data,2);
    Co=sum(Data,2);
    
    % Get peak threshold (minimum 2)
    ID=find(Co>=P_Th);
    
    % Get weighted adjacency from peaks (at least 2 active cells in the same frame)
    WAdjacency=(Data(ID,:)'*Data(ID,:)).*(1-eye(C));
    WMax=max(WAdjacency(:));
    
    % Get weight threshold (minimum 1)
    if (WMax>0)
        if W_Th_Auto
            N=100; alpha=0.05;
            W_Th = WeightTh(Data, N, alpha);
            if isempty(W_Th)
                Adjacency=zeros(C);
            else
                Adjacency=double(WAdjacency>=W_Th);
            end
        end
    else
        W_Th=-1;
        Adjacency=zeros(C);
    end
    % Get binary adjacency matrix 
end