% Plot Network Resilience

function [h Links Clocal]=Plot_Resilience(NetworkName,Adjacency,idxCardinal)
 
    % Resilience
    N=length(Adjacency);
    idxRand=randperm(N);

    % Remove neurons randomly
    for i=1:N-1
        [MatrixWithoutE]=DelElementMatrix(Adjacency,idxRand(1:i));
        E1(i)=efficiency_bin(double(MatrixWithoutE));
        Comp1(i)=max(get_components(MatrixWithoutE));
    end

    % Remove neurons by cardinality index
    for i=1:N-1
        [MatrixWithoutE]=DelElementMatrix(Adjacency,idxCardinal(1:i));
        E2(i)=efficiency_bin(MatrixWithoutE);
        Comp2(i)=max(get_components(MatrixWithoutE));
    end

    % Resilience
    Set_Figure(['Resilience (' NetworkName ')'],[0 0 600 600]);
    
    xl=[0 100];
    x=(1:N-1)/N*100;
    
    Set_Axes('Resilience-E',[0 0.5 1 0.5])
    plot(x,E1,'ok','linewidth',2,'markersize',10);hold on
    plot(x,E2,'or','linewidth',2,'markersize',10)
    legend({'by random','by cardinality'},'location','best')
    title('Resilience');xlabel('% neurons removed'); ylabel('Efficiency')
    set(gca,'box','off','xlim',xl)

    Set_Axes('Resilience-E',[0 0 1 0.5])
    plot(x,Comp1,'ok','linewidth',2,'markersize',10);hold on
    plot(x,Comp2,'or','linewidth',2,'markersize',10)
    title('Resilience');xlabel('% neurons removed'); ylabel('Components')
    set(gca,'box','off','xlim',xl)
    
    
end