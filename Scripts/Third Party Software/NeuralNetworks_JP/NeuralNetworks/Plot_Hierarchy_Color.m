% Plot C(k) vs k WITH colors


function Experiment=Plot_Hierarchy_Color(Experiment)
    
    NetworkName=Experiment.Data.DataName;
    Adjacency=Experiment.Adjacency.AdjacencyOnlyConnected;
    XYnodecolor=Experiment.Network.XYnodecolor;
    idxConnected=Experiment.Network.idxConnected;

    % Number of cells
    N=length(Adjacency);
    
    if ~N
        h=0;
        Links=0;
        Clocal=0;
        return;
    end

    Links=sum(Adjacency)';
    Clocal=clustering_coef_bu(Adjacency);
    l=Links;
    c=Clocal;
    
   %{
    % Plot local coefficient mean
    
    k=1;
    c=[]; error=[]; l=[];
    maxLinks=max(Links);
    for i=1:maxLinks
        idx=find(Links==i);
        if (~isempty(idx))
            c(k)=mean(Clocal(idx));
            error(k)=std(Clocal(idx))/sqrt(length(idx));
            l(k)=i;
            k=k+1;
        end        
    end
    k=length(l);
    %}
    
    % Delete nodes with only one link
    idx=find(l>1);
    l=l(idx);
    c=c(idx);
    
    N=length(l);
    
    h=Set_Figure([NetworkName ' - Hierarchical organization test'],[0 0 600 300]);
    
    x=l;
    y=c;

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
    %ylabel('C(k)')
    text(2,.2,{['\alpha=' num2str(slope,'%0.3f')] ['R^2=' num2str(R2,'%0.3f')]},...
        'color',[0 0 0],'FontSize',14,'FontWeight','bold')
    ylim([.1 2]); xlim([1 100])
    xlabel('k'); %ylabel('C(k)')
    
    axes(findobj('tag',['Axes - C(k) lineal ' NetworkName]))
    for i=1:N
        plot(l(i),c(i),'.k','MarkerSize',35)
        plot(l(i),c(i),'.','MarkerSize',23,'color',XYnodecolor(idxConnected(i),:))
    end
    ylim([0 1.1]); %xlim([1 30])
    xlabel('k'); ylabel('C(k)')
    title('Local Clustering Coeficient C(k)')

    axes(findobj('tag',['Axes - C(k) loglog ' NetworkName]))
    for i=1:N
        loglog(l(i),c(i),'.k','MarkerSize',35)
        loglog(l(i),c(i),'.','MarkerSize',23,'color',XYnodecolor(idxConnected(i),:))
    end


    % Plot errors
    %{
    for i=1:k
        x=[l(i) l(i)];
        y=[c(i)-error(i) c(i)+error(i)];
        axes(findobj('tag',['Axes - C(k) lineal' NetworkName]))
        plot(x,y,'color',[1 .8 .8],'linewidth',2,'markersize',10)
        axes(findobj('tag',['Axes - C(k) loglog' NetworkName]))
        loglog(x,y,'color',[1 .8 .8],'linewidth',2,'markersize',10)
    end
    %}
    
    Experiment.Hierarchy.Slope=slope;
    Experiment.Hierarchy.R2=R2;
    Experiment.Hierarchy.XYreal=[x y];
    Experiment.Hierarchy.XYfit=[xfit; yfit];
end