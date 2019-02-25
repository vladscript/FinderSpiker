% Plot network
% 
% Plot in circle
% 
% PlotNetwork(G, Fig)
%
% Inputs:
% G = binary undirected connection matrix
% Fig = number of figure
%
% ..:: by Jesús E. Pérez-Ortega ::.. March-2013

function PlotNetwork_JP(G,Fig)
    N=length(G);
    
    % Get XY coordinates of circle
    XY=[cos(2*pi*[1:N]'/N) sin(2*pi*[1:N]'/N)];
    
    % Color and width
    color_edge=[0 0 0];
    color_node=[.9 .9 .9];
    BGColor=[.8 .8 .8];
    width=2;
    
    % Plot edges
    if Fig
        figure(Fig)
    else
        figure('name','Network','position',[0 0 600 600]);
    end
    clf
    for a=1:N-1
        for b=a+1:N
            if G(a,b)
                plot(XY([a b],1),XY([a b],2),'color',color_edge,'LineWidth',width);
                hold on
            end
        end
    end
    
    % Plot neurons
    for i=1:N
        plot(XY(i,1),XY(i,2),'.k','MarkerSize',35)
        plot(XY(i,1),XY(i,2),'.','color',color_node,'MarkerSize',23)
        text(XY(i,1)*1.1,XY(i,2)*1.1,num2str(i),'HorizontalAlignment','Center')
    end
    set(gca,'view',[0 -90])
    set(gca,'xlim',[-1.2 1.2],'ylim',[-1.2 1.2])
    set(gca,'XTickLabel',[],'YTickLabel',[],'color',BGColor,'XColor',BGColor,'YColor',BGColor)
end