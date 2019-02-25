% Verifies if Network satisfies the (ultra) Small-Wrold properties
%
% Jesús Pérez Apr-2014

function [IsUltra IsSmallWorld Properties] = SmallWorld_JP(Experiment)
    
    NetworkName=Experiment.Data.DataName;
    Adjacency=Experiment.Adjacency.AdjacencyOnlyConnected;


    Set_Figure(['SmallWorld - ' NetworkName],[0 0 600 300]);

    % Total nodes and total edges
    N = length(Adjacency);        
    K = sum(sum(Adjacency))/2;    
    
    %% Comprobar que la red es de escala libre
    [x y xfit yfit slope R2] = FitPowerLaw(Adjacency);
    
    % Plot Degree distribution
    Set_Axes(['DegreeDistribution - ' NetworkName],[0 0 .5 1 ])
    ColorLight=[0.9 0 0.7];
    ColorDark=[0.5 0 0.3];
    Green=[0 .6 0];
    Red=[.8 0 0];
    plot(x,y,'o','color',ColorLight); hold on;
    plot(xfit,yfit,'--','color',ColorDark)
    text(max(x)/2,max(y)*0.8,['\alpha=' num2str(-slope,'%0.1f')],...
        'FontSize',14,'FontWeight','bold')
    text(max(x)/2,max(y)*0.6,['R^2=' num2str(R2,'%0.3f')],...
        'FontSize',14,'FontWeight','bold')
    title('Degree Distribution'); xlabel('# neuron'); ylabel('Count of links')    
    
    xtext=xlim;
    ytext=ylim;
    if (R2 > 0.5)
        IsUltra=true;
        text(xtext(1),ytext(1),'Ultra','color',Green,'FontSize',24,...
            'VerticalAlignment','bottom','FontWeight','bold')
    else
        IsUltra=false;
        text(xtext(1),ytext(1),'X','color',Red,'FontSize',24,...
            'VerticalAlignment','bottom','FontWeight','bold')
    end
    
    %% Comprobar que la red tiene propiedades de mundo pequeño

    % Red real
    % Cálculo de la longitud característica
    D=distance_bin(Adjacency);              
    [L E]=charpath(D);                      

    % Cálculo del coeficiente de agrupamiento
    Clocal=clustering_coef_bu(Adjacency);   
    C=mean(Clocal);                         
    
    % Red regular
    Reg=RegNet_JP(N,K);
    % Cálculo de la longitud característica
    Dreg=distance_bin(Reg);                     
    [Lreg Ereg]=charpath(Dreg);                 
    % Cálculo del coeficiente de agrupamiento
    Clocalreg=clustering_coef_bu(Reg);          
    Creg=mean(Clocalreg);                       

    % Red aleatoria
    Rand=makerandCIJ_und(N,K);
    % Cálculo de la longitud característica
    Drand=distance_bin(Rand);                   
    [Lrand Erand]=charpath(Drand);              
    % Cálculo del coeficiente de agrupamiento
    Clocalrand=clustering_coef_bu(Rand);        
    Crand=mean(Clocalrand);
    % Medida de small-world
    Omega=Lrand/L-C/Creg;
    
    minC=min([Creg Crand]);
    maxC=max([Creg Crand]);
    
    minL=min([Lreg Lrand]);
    maxL=max([Lreg Lrand]);
    
    minE=min([Ereg Erand]);
    maxE=max([Ereg Erand]);
    
    % Plot figure 
    Set_Axes(['Network Properties - ' NetworkName],[.5 0 .5 1 ])
    
    normC=([Creg C Crand]-minC)/(maxC-minC);
    normL=([Lreg L Lrand]-minL)/(maxL-minL);
    normE=([Ereg E Erand]-minE)/(maxE-minE);
    
    plot(normC,'-ob'); hold on;
    plot(normL,'-or')
    plot(normE,'--or'); hold off;
    
    set(gca,'xtick',1:3,'xticklabel',{'Regular','Real','Random'})
    xlim([.5 3.5])
    title(['Properties: L=' num2str(L,'%1.2f') ' E=' num2str(E,'%1.2f')...
        ' C=' num2str(C,'%1.2f') ' Omega=' num2str(Omega,'%1.2f')])
    legend('C', 'L', 'E')
    
    xtext=xlim;
    ytext=ylim;
    
    normC=normC(2);
    normL=normL(2);
    normE=normE(2);
    
    if (normC>0.9 && normL<0.1)
        IsSmallWorld=true;
        text(xtext(1),ytext(1),'Small World','color',Green,'VerticalAlignment',...
            'bottom','FontSize',24,'FontWeight','bold')
    else
        IsSmallWorld=false;
        text(xtext(1),ytext(1),'X','color',Red,'VerticalAlignment',...
            'bottom','FontSize',24,'FontWeight','bold')
    end
    
    Properties.N=N;
    Properties.K=K;
    
    Properties.slope=-slope;
    Properties.R2_slope=R2;
    
    Properties.C=C;
    Properties.L=L;
    Properties.E=E;
    
    Properties.Crand=Crand;
    Properties.Lrand=Lrand;
    Properties.Erand=Erand;
    
    Properties.Cnorm=normC;
    Properties.Lnorm=normL;
    Properties.Enorm=normE;
    
    Properties.Omega=Omega;
    
end