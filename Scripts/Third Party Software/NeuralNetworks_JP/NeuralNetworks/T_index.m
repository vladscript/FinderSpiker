%% T_index  T index  
    %%which is the clustering index, and obtains the value of it.

%A is the input variable that recibes any Matrix 

%T is the variable that contains the value of the cluster index

% Jesús E. Pérez-Ortega (modified 08-12)


function [T]=T_index(A)

%start Values
totalT=0; % the total actual value of triangles in the newtwork  
Triang=0; % the total possible triangles in the network

%cicle that calculates the index T
x=length(A);
for i=1:x
    for j=1:x-1
        for m=j+1:x       % works with the row i
            if(i~=m)       %   discards the matrix diagonal
                if(A(i,j)&& A(i,m)) % to find a possible triangle
                    totalT=totalT+1;
                    if (A(j,m)) % confirm a triangle
                       Triang=Triang+1;
                    end
                end
            end
        end
    end
end
T=Triang/totalT;