%% Geodesic distance
% JP Oct-12

function [Geo GeoSingleMean GeoAllMean Cuts]=Geodesic_JP(Adjacent)
    Adjacent=double(Adjacent);
    N=length(Adjacent);
    
    % Geodesic mean
    Geo=MinPath(Adjacent);
    
    % Geodesic mean of single node
    for i=1:N
        single=Geo(:,i);
        GeoSingleMean(i)=mean(single(single~=inf&single~=0));
    end
    
    % Geodesic mean of the network
    TriGeo=triu(Geo,1);
    TriGeoClean=TriGeo(TriGeo~=inf&TriGeo~=0);
    GeoAllMean=mean(TriGeoClean);
    
    % Cuts of the network
    Cuts=length(find(TriGeo==inf))/N^2*2;
end

function MinPathMatrix=MinPath(Adjacent)
% mindist_wpau  Distance
% returns the minimum distance of the Network  for each node 
%A is the input variable that recibes any Matrix
%A (out) is the output matrix that has all the distances to each node
%B is the matrix that contains the number of iteration (or node) where it
%  found the minimum distance

% JP modified 30-sep-12

     N=size(Adjacent,1);
         
     Adjacent(Adjacent==0)=inf;
     for i=1:N
         Adjacent(i,i)=0;
     end
     MinPathMatrix=Adjacent;
     
     for k=1:N
        for i=1:N-1
            for j=i+1:N
               if(k~=i)
                   if (k~=j)
                       c=MinPathMatrix(k,i)+MinPathMatrix(k,j); % the sum of the positions in the weigth matrix
                       if(c<MinPathMatrix(i,j)) %the value of c has to be less than the value of the position i,j
                           MinPathMatrix(i,j)=c; %then it sets the value
                           MinPathMatrix(j,i)=c;
                       end
                   end
               end
            end
        end
     end
end

%{
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
%}