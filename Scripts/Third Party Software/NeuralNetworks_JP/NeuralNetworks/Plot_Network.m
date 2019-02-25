%% Plot Network
function h=Plot_Network(NetworkName,Adjacency,XY,XY_Colors)
    
    if ~exist('XY_Colors')
        XY_Colors=ones(length(Adjacency),3);
    end
        
    h=Set_Figure(NetworkName,[0 0 600 600]);
    Set_Axes(['Axes' NetworkName],[0 0 1 1])
    Axes_White();
    Plot_EdgesAndNodes(Adjacency,XY,XY_Colors);
    Axes_XYFit(XY);
    Axes_Margin()
    title(NetworkName)
end