%% Set Figure
function h=Set_Figure(TitleName, Position)
    h=findobj('name',TitleName);
    if isempty(h)
        h=figure('name',TitleName,'numbertitle','off','position',Position);
    else
        figure(h);
    end
    clf
end
