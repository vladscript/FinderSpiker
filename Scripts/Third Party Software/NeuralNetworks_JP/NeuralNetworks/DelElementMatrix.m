function [MatrixWithoutE]=DelElementMatrix(Matrix,E)
% JP 30-sep-12
    N=length(Matrix);
    Es=sort(E);
    k=length(E);
    for i=1:k
        E=Es(i)-i+1;
        if E==1
            MatrixWithoutE=Matrix(2:N,2:N);
        elseif E==N
            MatrixWithoutE=Matrix(1:N-1,1:N-1);
        else
            MatrixWithoutE=[Matrix(1:E-1,1:E-1) Matrix(1:E-1,E+1:N);...
                Matrix(E+1:N,1:E-1) Matrix(E+1:N,E+1:N)];
        end
        N=length(MatrixWithoutE);
        Matrix=MatrixWithoutE;
    end
    MatrixWithoutE=double(MatrixWithoutE);
end
