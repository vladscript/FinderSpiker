% Generate Hierarchical Network
%
% This function generates binary undirected network with hierarchical
% connections.
%
% H = HNet(N,K,C)
%
% Inputs:
% N = number of nodes
% K = number of links
% C = clustering coefficient
%
% Outputs:
% H = binary undirected adjacent matrix of hierarchical network
%
% ..:: by Jesús E. Pérez-Ortega ::.. Sep-2014

function H = HNet(N,K,C)
    
    % Enlaces totales
    %K=round(rho*N^2-N);

    Final=zeros(N);

    m=3;
    rep=round(m*C);

    %% Paso 1
    pattern=ones(m)-eye(m);
    Final(1:m,1:m)=pattern;

    % actualiza probabilidad
    nLinks=sum(Final(:));
    Pnode=sum(Final)/nLinks;

    %% Paso 2
    % Jerarquía 1
    for i=(m+1):m:(N-m+1)

        % agrega m=3 nodos
        id=i:(i+m-1);
        Final(id,id)=pattern;


        for j=1:rep
            for k=2:(m)
                if (rand<C)
                    % conecta m-1 nodos
                    new=id(k);

                    % Selecciona el 1er nodo a enlazar de acuerdo a la probabilidad de unión
                    [sortPnode idx]=sort(Pnode,'descend');

                    n=idx(1+sum(rand>cumsum(sortPnode)));
                    %
                    while(Final(new,n)) 
                        n=idx(1+sum(rand>cumsum(sortPnode)));
                    end
                    %}
                    Final(new,n)=1;
                    Final(n,new)=1;
                end

            end

            % actualiza probabilidad
            nLinks=sum(Final(:))/2;
            Pnode=sum(Final)/nLinks;
        end

    end
    %% Paso 3
    % Jerarquía 2
    for j=2:floor(log(N)/log(m))
        size_p=m^j;
        pattern=Final(1:size_p,1:size_p);
        for i=(size_p+1):size_p:(N-size_p-1)
            if (rand<C)
                id=i:(i+size_p-1);
                Final(id,id)=pattern;
            end
        end
    end

    % actualiza probabilidad
    nLinks=sum(Final(:))/2;
    Pnode=sum(Final)/nLinks;

    %% Paso 4
    notConnected=find(sum(Final)==0);
    if (notConnected)
        num=length(notConnected);
        for i=1:num
            n1=notConnected(i);
            [sortPnode idx]=sort(Pnode,'descend');
            n2=idx(1+sum(rand>cumsum(sortPnode)));
            Final(n1,n2)=1;
            Final(n2,n1)=1;
        end
    end

    % actualiza probabilidad
    nLinks=sum(Final(:))/2;
    Pnode=sum(Final)/nLinks;

    %% Paso 5

    links=sum(Final(:))/2;
    fit=K-links;

    while (fit>0)
        n1=1+round(rand*(N-1));

        [sortPnode idx]=sort(Pnode,'descend');
        n2=idx(1+sum(rand>cumsum(sortPnode)));

        Final(n1,n2)=1;
        Final(n2,n1)=1;

        links=sum(Final(:))/2;
        fit=K-links;
    end
    %
    mix=randperm(N);
    Final=Final(mix,mix);
    %}
    H=Final;
%{
    %% 6
    K=sum(Final(:)/2);
    k=K/N;
    figure(1); clf;
    imagesc(Final)
    %
    Adjacency=Final;
    N = length(Adjacency);        
    K = sum(sum(Adjacency))/2;
    k=K/N;
    %

    Experiment.Data.DataName='Model H';
    Experiment.Adjacency.AdjacencyOnlyConnected=Adjacency;
    SmallWorld_JP2(Experiment);


    Adjacency=Final;

    
    C=clustering_coef_bu(Adjacency);
    links=sum(Adjacency);

    x=links;
    y=C;

    figure(200);clf;
    subplot(2,1,1)
    plot(x,y,'or')
    subplot(2,1,2)
    loglog(x,y,'or');hold on

    N=length(Adjacency);
    [slope, intercept, R2] = logfit(x,y,'loglog');
    xfit=min(x):(max(x)-min(x))/100:max(x);
    yfit=(10^intercept)*xfit.^(slope);
    loglog(xfit,yfit,'--','linewidth',2.5)
    %}

end