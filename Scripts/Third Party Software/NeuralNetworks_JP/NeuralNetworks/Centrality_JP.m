%% AveDstWN_PR Average Distance Without a Node in the Network (Pau Ramírez)
%returns 2 graphs, a Histogram  that contains all the changes in path of
%the newtwork without a node(a subgraph)
%and a Simple Plot with a mean value
%Also this function converts the recipt matrix to a double matrix so it can
%be useful to the own function, if it is not, the function doesn't work .

%Adjacent is the Adyancecy Matrix of the Network
% TrimsInNet is an array that contains the following ChangesInMinPaths/tot, wich is a mean
% value

%ChangesInMinPaths is an array that contains the total changes in paths in the network caused for the disapperance of
%a node
%tot is an array that contains the total valid paths in the subgraph
%



% JP modified 30-sep-12
function [ChangesInMinPaths TrimsInNet]=Centrality_JP(Adjacent)
    Adjacent=double(Adjacent);
    MinReal=MinPath(Adjacent);
    
    N=length(MinReal);
    
    ChangesInMinPaths=zeros(1,N);
    TrimsInNet=zeros(1,N);
    for i=1:N
        [AjacentComp]=DelElementMatrix(Adjacent,i);
        [MinRealComp]=DelElementMatrix(MinReal,i);
        MinDel=MinPath(AjacentComp);
        ChangesInMinPaths(i)=length(find(MinRealComp~=MinDel));
        TrimsRealComp=length(find(MinRealComp==inf))/2;
        TrimsInNet(i)=length(find(MinDel==inf))/2-TrimsRealComp;
    end
end

function [MatrixWithoutE]=DelElementMatrix(Matrix,E)
% JP 30-sep-12
    N=length(Matrix);
    if E==1
        MatrixWithoutE=Matrix(2:N,2:N);
    elseif E==N
        MatrixWithoutE=Matrix(1:N-1,1:N-1);
    else
        MatrixWithoutE=[Matrix(1:E-1,1:E-1) Matrix(1:E-1,E+1:N);...
            Matrix(E+1:N,1:E-1) Matrix(E+1:N,E+1:N)];
    end
 end