% Randomize regular network
%
% RandNet = RandRegNet_JP(Reg, P)
%
% Inputs:
% Reg = regular network as N x N matrix (N = #nodes)
% P = probability of rewiring (0 to 1)
%
% Outputs:
% RandNet = randomized network wit P probability of rewiring
%
% ..:: by Jesús E. Pérez-Ortega ::.. March-2013
% Modified randi instead of rand

function RandNet = RandRegNet_JP(Reg, P)
    if P
        N=length(Reg);
        for i=1:N
            idx=find(Reg(i,:));
            nj=length(idx);
            for j=1:nj
                % keep=rand(round(1/P))-1;
                keep=randi(round(1/P),1)-1;
                if ~keep
                    Reg=Reg+eye(N);
                    no_edge=find(Reg(i,:)==0);
                    new=randi(length(no_edge));
                    new_edge=no_edge(new);
                    while new_edge==i
                        new=randi(length(no_edge));
                        new_edge=no_edge(new);
                    end
                    Reg(i,new_edge)=1;
                    Reg(new_edge,i)=1;
                    Reg(i,idx(j))=0;
                    Reg(idx(j),i)=0;
                    Reg=Reg-eye(N);
                end
            end
        end
    end
    RandNet=Reg;
end