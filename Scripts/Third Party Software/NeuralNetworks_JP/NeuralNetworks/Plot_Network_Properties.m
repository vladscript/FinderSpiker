% Plot Network with Properties

function [h h2]=Plot_Network_Properties(Experiment)
    
    
    NetworkName=Experiment.Data.DataName;
    Adjacency=Experiment.Adjacency.AdjacencyOnlyConnected;
    IdxCardinal=Experiment.Ranks.Cardinality;
    
    % Number of cells
    N=length(Adjacency);
    if ~N
        h=0;
        h2=0;
        return;
    end

    %
    % Generate XY circular coordinates
    XY=Get_XYCircular(N);
    
    % Plot Network
    %h=Plot_Network(NetworkName,Adjacency,XY)
    h=Plot_ReorderedNetwork(NetworkName, Adjacency);
    
    % Links density
    K=sum(Adjacency(:))/2;
    Rho=K/(N*(N-1));

    % Characteristic path lenght and Efficiency
    D=distance_bin(Adjacency);
    [L E]=charpath(D);
    
    % Clustering coefficient
    c=clustering_coef_bu(Adjacency);
    C=mean(c);

    % Assotativity
    A=assortativity_newman(Adjacency);
    
    title([NetworkName ' - N=' num2str(N) ' K=' num2str(K)...
        'rho=' num2str(Rho,'%1.2f') '  L=' num2str(L,'%1.2f')...
        '  C=' num2str(N,'%1.2f') 'C2=' num2str(C,'%1.2f') '  E=' num2str(E,'%1.2f')...
        '  A=' num2str(A,'%1.2f') ])
    
    % Degree distribution
    %h2=Plot_DegreeDistribution(NetworkName,Adjacency);
    
    % Plot Resilience
    Plot_Resilience(NetworkName,Adjacency,IdxCardinal);
    
    %}
end

