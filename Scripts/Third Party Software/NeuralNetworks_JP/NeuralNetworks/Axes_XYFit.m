%% Axes XY Fit
function Axes_XYFit(XY)
    Xmin=min(XY(:,1));
    Xmax=max(XY(:,1));
    Ymin=min(XY(:,2));
    Ymax=max(XY(:,2));
    xlim([Xmin Xmax])
    ylim([Ymin Ymax])
end