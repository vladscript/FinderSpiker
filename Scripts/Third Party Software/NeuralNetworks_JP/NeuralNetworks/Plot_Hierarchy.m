% Plot C(k) vs k

function [h Links Clocal]=Plot_Hierarchy(Experiment)
    
    NetworkName=Experiment.Data.DataName;
    Adjacency=Experiment.Adjacency.AdjacencyOnlyConnected;

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
    
   %
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
        
    % delete zeros
    l=l(find(c>0));
    c=c(find(c>0));
    
    
    h=Set_Figure([NetworkName ' - Hierarchical organization test (mean+SEM)'],[0 0 600 600]);
    
    x=l;
    y=c;

    [slope, intercept, R2] = logfit(x,y,'loglog');
    xfit=min(x):(max(x)-min(x))/100:max(x);
    yfit=(10^intercept)*xfit.^(slope);

    Set_Axes(['Axes - C(k) lineal' NetworkName],[0 0.5 1 0.5]); hold on
    plot(xfit,yfit,'--','linewidth',2,'markersize',10,'color',[1 .8 .8])

    Set_Axes(['Axes - C(k) loglog' NetworkName],[0 0 1 0.5])
    loglog(xfit,yfit,'--','linewidth',2,'markersize',10,'color',[1 .8 .8]); hold on
    set(gca,'tag',['Axes - C(k) loglog' NetworkName])
    xlabel('k')
    ylabel('C(k)')
    text(2,.2,['slope=' num2str(slope,'%0.3f') '      R^2=' num2str(R2,'%0.3f')],...
        'color',[.8 .3 .3],'FontSize',14,'FontWeight','bold')
    ylim([.1 1]); xlim([1 100])
    xlabel('k'); ylabel('C(k)')
    
    axes(findobj('tag',['Axes - C(k) lineal' NetworkName]))
    plot(l,c,'or','linewidth',2,'markersize',10)
    ylim([0 1]); %xlim([1 30])
    xlabel('k'); ylabel('C(k)')
    title('Local Clustering Coeficient C(k)')

    axes(findobj('tag',['Axes - C(k) loglog' NetworkName]))
    loglog(l,c,'or','linewidth',2,'markersize',10)

    % Plot errors
    %
    for i=1:k
        x=[l(i) l(i)];
        y=[c(i)-error(i) c(i)+error(i)];
        axes(findobj('tag',['Axes - C(k) lineal' NetworkName]))
        plot(x,y,'color',[1 .8 .8],'linewidth',2,'markersize',10)
        axes(findobj('tag',['Axes - C(k) loglog' NetworkName]))
        loglog(x,y,'color',[1 .8 .8],'linewidth',2,'markersize',10)
    end
    %}

end