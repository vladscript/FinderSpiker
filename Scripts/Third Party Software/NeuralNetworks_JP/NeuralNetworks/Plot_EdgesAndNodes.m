%% Plot Edges And Nodes
function Plot_EdgesAndNodes(Adjacency,XY,XY_Colors)
    
    Adjacency=double(Adjacency>0);
    
    % Plot edges
    width=1;
    C=length(Adjacency);
    color=[0 0 0];
    
    hold on
    for a=1:C-1
        for b=a+1:C
            if (Adjacency(a,b))
                plot(XY([a b],1),XY([a b],2),'color',color,'LineWidth',width);
            end
        end
    end

    % Plot nodes
    color=[1 1 1];
    for i=1:C
        plot(XY(i,1),XY(i,2),'.k','MarkerSize',35)
        plot(XY(i,1),XY(i,2),'.','color',XY_Colors(i,:),'MarkerSize',23)
        text(XY(i,1)*1.1,XY(i,2)*1.1,num2str(i),'HorizontalAlignment','Center')
    end
    
end