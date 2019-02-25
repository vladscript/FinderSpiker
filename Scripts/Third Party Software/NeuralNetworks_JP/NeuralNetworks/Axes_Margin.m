%% Axes Margin
function Axes_Margin()
    Margin=0.05;
    X=get(gca,'xlim');
    Y=get(gca,'ylim');
    Xoff=(X(2)-X(1))*Margin;
    Yoff=(Y(2)-Y(1))*Margin;
    xlim([X(1)-Xoff X(2)+Xoff]);
    ylim([Y(1)-Yoff Y(2)+Yoff]);
end