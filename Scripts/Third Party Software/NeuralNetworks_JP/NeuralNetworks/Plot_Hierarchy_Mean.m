% Plot Network with Properties

function Experiment=Plot_Hierarchy_Mean(Experiment)

    NetworkName=Experiment.Data.DataName;
    Adjacency=Experiment.Adjacency.AdjacencyOnlyConnected;
    
    
    Set_Figure([NetworkName ' - Hierarchical organization test'],[0 0 600 300]);
    Set_Axes(['Axes - Hierarchical' NetworkName],[0 0 1 1])
    
    % Number of cells
    Links=sum(Adjacency)';
    Clocal=clustering_coef_bu(Adjacency);
        
    % Delete nodes with only one link
    idx=find(Links>1);
    Links=Links(idx);
    Clocal=Clocal(idx);
    
    %
    maxLinks=max(Links);
    
    k=1;
    for i=1:maxLinks
        idx=find(Links==i);
        if (~isempty(idx))
            c(k)=mean(Clocal(idx));
            error(k)=std(Clocal(idx))/sqrt(length(idx));
            l(k)=i;
            k=k+1;
        end        
    end
    Links=l;
    Clocal=c;
    %}
    
    Set_Figure([NetworkName ' - Hierarchical organization test'],[0 0 600 300]);
    N=length(Links);
    x=Links;
    y=Clocal;

    [slope, intercept, R2] = logfit(x,y,'loglog');
    xfit=min(x):(max(x)-min(x))/100:max(x);
    yfit=(10^intercept)*xfit.^(slope);

    Set_Axes(['Axes - C(k) lineal ' NetworkName],[0 0 0.5 1]); hold on
    plot(xfit,yfit,'--k','linewidth',2,'markersize',10)
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    set(gca,'box','off')

    Set_Axes(['Axes - C(k) loglog ' NetworkName],[0.5 0 0.5 1])
    loglog(xfit,yfit,'--k','linewidth',2,'markersize',10); hold on
    set(gca,'FontSize',14,'LineWidth',2,'TickLength',[0.02 0.02])
    set(gca,'box','off')
    set(gca,'tag',['Axes - C(k) loglog ' NetworkName])
    xlabel('k')
    ylabel('C(k)')
    text(2,.2,{['\alpha=' num2str(slope,'%0.3f')] ['R^2=' num2str(R2,'%0.3f')]},...
        'color',[0 0 0],'FontSize',14,'FontWeight','bold')
    ylim([.1 2]); xlim([1 100])
    xlabel('k'); ylabel('C(k)')
    
    % Plot lineal
    axes(findobj('tag',['Axes - C(k) lineal ' NetworkName]))
    
    % Plot errors
    for i=1:N
        x=[l(i) l(i)];
        y=[c(i)-error(i) c(i)+error(i)];
        plot(x,y,'k','linewidth',2)
    end
    for i=1:N
        plot(Links(i),Clocal(i),'.k','MarkerSize',35)
        plot(Links(i),Clocal(i),'.','MarkerSize',23,'color',[0 0 .9])
    end
    ylim([0 1.1]); %xlim([1 30])
    xlabel('k'); ylabel('C(k)')
    title('Local Clustering Coeficient C(k)')

    % Plot loglog
    axes(findobj('tag',['Axes - C(k) loglog ' NetworkName]))
    for i=1:N
        x=[l(i) l(i)];
        y=[c(i)-error(i) c(i)+error(i)];
        loglog(x,y,'k','linewidth',2)
    end
    for i=1:N
        loglog(Links(i),Clocal(i),'.k','MarkerSize',35)
        loglog(Links(i),Clocal(i),'.','MarkerSize',23,'color',[0 0 .9])
    end
    
    Experiment.Hierarchy.Slope=slope;
    Experiment.Hierarchy.R2=R2;
    Experiment.Hierarchy.XYreal=[x y];
    Experiment.Hierarchy.XYfit=[xfit; yfit];
end