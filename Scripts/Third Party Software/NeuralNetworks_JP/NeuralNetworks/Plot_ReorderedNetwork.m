% Plot Reordered Network
function h=Plot_ReorderedNetwork(DataName, Adjacency)
    
    % Get Connected
    [Connected Disconnected NumCon NumDis]=Get_Connection(Adjacency);
    XYLine=Get_XYLinear(NumDis);
    XYCirc=Get_XYCircular(NumCon);
    
    % Reorder
    if max(Adjacency(:))
        [IndicesRe Adj Cost Costs Epoch]=ReorderMAT_JP(Adjacency(Connected, Connected),10000,'circ');
        XYCircReorder(IndicesRe,:)=XYCirc;
        XY([Connected Disconnected],:)=[XYCircReorder;XYLine];
    else
        XY=XYLine;
    end
    
    % Plots
    h=Plot_Network('Reordered Network',Adjacency,XY);
    title(['Reordered Network - ' DataName ' (Epochs: ' num2str(Epoch) ' - Cost: ' num2str(round(Cost)) ')'])

end