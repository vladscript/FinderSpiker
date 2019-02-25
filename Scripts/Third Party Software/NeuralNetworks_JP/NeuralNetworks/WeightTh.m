% Weight threshold
% 
% Montecarlo for determining weight threshold in weighted matrix.
% 
% Threshold is determined when the edges number is significantly higher in
% original net versus its randomized version.
%
% Adjacency matrix is determined by synchrony, weight of connection means
% the times which are synchronized the nodes linked.
%
% W_Th = WeightTh(Data, N, alpha)
% 
% ..:: by Jesús E. Pérez-Ortega ::.. Jun-2013 

function W_Th = WeightTh(Data, N, alpha)

    % Data size
    [F C]=size(Data);

    % Adjacent from original net
    A=Data'*Data.*(1-eye(C));

    % Number of edges at each threshold
    edges_th_A=[];
    for Ths=1:F
        edges_th_A(Ths)=sum(sum(double(A>Ths)))/2;
    end

    % Higher values in original net
    higher=zeros(1,F);

    for j=1:N
        % Randomized version
        DataRnd=[];
        for i = 1:C
            DataRnd(:,i) = Data(randperm(F),i);  % Shuffle spikes for every cell
        end

        % Adjacent from randomized version
        ARnd=DataRnd'*DataRnd.*(1-eye(C));

        % Number of edges at each threshold
        edges_th_ARnd=[];
        for Ths=1:F
            edges_th_ARnd(Ths)=sum(sum(double(ARnd>Ths)))/2;
        end

        % Times where edges from original net are higher than randomized one
        higher_j=double((edges_th_A-edges_th_ARnd)>0);
        higher=higher+higher_j;
    end

    % Weight threshold
    [value W_Th]=find(higher>((1-alpha)*N),1,'first');
end