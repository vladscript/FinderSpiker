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