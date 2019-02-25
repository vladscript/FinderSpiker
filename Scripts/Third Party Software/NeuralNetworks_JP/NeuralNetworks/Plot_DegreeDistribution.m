% Plot degree distribution

function h=Plot_DegreeDistribution(NetworkName, Adjacency)
    
    C=length(Adjacency);
    Links=sum(Adjacency,1);

    if ~Links
        return;
    end
    
    % Colors
    BGColor=[1 1 1];
    ColorLight=[0.9 0 0.7];
    ColorDark=[0.5 0 0.3];
    ColorTextLines=[0 0 0];
    
    % Plot Links
    h=Set_Figure(['Links - ' NetworkName],[0 0 600 300]);
    Set_Axes(['Links ' NetworkName],[0 0 1 1])

    plot(Links,'-k')
    set(gca,'xlim',[0 C],'ylim',[0 max(Links)+1])
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    title(NetworkName,'color',ColorTextLines)
    xlabel('# neuron')
    ylabel('Count of links')
    set(gca,'Tag',['LinksAxes' NetworkName])
    
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
    plot(x,y,'.','color',ColorLight,'MarkerSize',30)
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    %[y x]=hist(Links,0:max_Links);
    %bar(x(2:max_Links+1),y(2:max_Links+1),'facecolor',ColorLight,'edgecolor',ColorDark)
    hold on
    set(gca,'color',BGColor,'box','off','xcolor',ColorTextLines,'ycolor',ColorTextLines)
    %title('Links histogram','color',ColorTextLines)
    xlabel('Links')
    ylabel('Count of neurons')

    % Plot power law fit (only neurons with at least one link)
    [slope, intercept, R2] = logfit(x,y,'loglog');
    xfit=min(x):(max(x)-min(x))/100:max(x);
    yfit=(10^intercept)*xfit.^(slope);
    plot(xfit,yfit,'--','color',ColorDark,'linewidth',2.5)
    % Write the function fitted
    text(max(x)/2,max(y)*0.8,['\alpha=' num2str(-slope,'%0.1f')],'color',ColorDark,'FontSize',14,'FontWeight','bold')
    text(max(x)/2,max(y)*0.6,['R^2=' num2str(R2,'%0.3f')],'color',ColorDark,'FontSize',14,'FontWeight','bold')

    % Loglog activity histogram
    Set_Axes(['HistLinksPL ' NetworkName],[0.5 0 0.5 1])
    loglog(x,y,'.','color',ColorLight,'MarkerSize',30)
    hold on
    loglog(xfit,yfit,'--','color',ColorDark,'linewidth',2.5)
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    set(gca,'box','off')
end