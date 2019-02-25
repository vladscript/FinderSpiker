% Plot degree distribution WITH color

function Experiment=Plot_DegreeDistribution_Color(Experiment)

    NetworkName=Experiment.Data.DataName;
    Adjacency=Experiment.Adjacency.AdjacencyOnlyConnected;
    XYnodecolor=Experiment.Network.XYnodecolor;
    idxConnected=Experiment.Network.idxConnected;
    
    C=length(Adjacency);
    Links=sum(Adjacency,1);

    if ~Links
        return;
    end
    
    % Colors
    BGColor=[1 1 1];
    ColorLight=[1 0.9 0.98];
    ColorDark=[0.5 0 0.3];
    ColorTextLines=[0 0 0];
    
    Set_Figure(['Distribution - ' NetworkName],[0 0 600 300]);

    % Plot links histogram
    Set_Axes(['Links ' NetworkName],[0 0 0.5 1])
    
    max_Links=max(Links);
    if ~max_Links
        return;
    end
    
    N=length(Links);
    nbins=round(log2(N)+1);
    [y x]=hist(Links,nbins);%This is a probe ,0:max_Links);
    y=y/N;
    bar(x,y,'facecolor',ColorLight,'edgecolor',ColorDark); hold on

    % Plot power law fit (only neurons with at least one link)
    [slope, intercept, R2] = logfit(x,y,'loglog');
    xfit=min(x):(max(x)-min(x))/100:max(x);
    yfit=(10^intercept)*xfit.^(slope);
    plot(xfit,yfit,'--','color',ColorDark,'linewidth',2.5)
    % Write the function fitted
    

    plot(x,y,'.','color',ColorDark,'MarkerSize',35)
    plot(x,y,'.','color',ColorLight,'MarkerSize',23)
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    title('Degree distribution P(k)','color',ColorTextLines)
    xlabel('k'); ylabel('P(k)')
    
    % Loglog activity histogram
    Set_Axes(['HistLinksPL ' NetworkName],[0.5 0 0.5 1])
    loglog(xfit,yfit,'--','color',ColorDark,'linewidth',2.5); hold on
    loglog(x,y,'.','color',ColorDark,'MarkerSize',35)
    loglog(x,y,'.','color',ColorLight,'MarkerSize',23)
    text(2,.1,{['\gamma=' num2str(-slope,'%0.1f')] ['R^2=' num2str(R2,'%0.3f')]},...
        'color',[0 0 0],'FontSize',14,'FontWeight','bold')
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    set(gca,'box','off')
    xlabel('k');
    
    
    Experiment.DegreeDistribution.Slope=slope;
    Experiment.DegreeDistribution.R2=R2;
    Experiment.DegreeDistribution.XYreal=[x; y];
    Experiment.DegreeDistribution.XYfit=[xfit; yfit];
end