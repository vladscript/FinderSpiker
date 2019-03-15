if exist('RoA_POS')
    aremerged=true;    % Are there already colocated cells
else
    aremerged=false;    % Are there already colocated cells
end

plot_cdf_cell(RoA_ALL,Names_Conditions);
xlabel('RoA')
plot_cdf_cell(ISIs_ALL,Names_Conditions);
xlabel('ISI')
plot_cdf_cell(TranLengths_ALL,Names_Conditions);
xlabel('ED')
if aremerged
    % POSITIVE
    plot_cdf_cell(RoA_POS,Names_Conditions);
    xlabel('RoA'); set(gcf,'Name','POSITIVE')
    plot_cdf_cell(ISIs_POS,Names_Conditions);
    xlabel('ISI'); set(gcf,'Name','POSITIVE')
    plot_cdf_cell(TranLengths_POS,Names_Conditions);
    xlabel('ED'); set(gcf,'Name','POSITIVE')
    % NEGATIVE
    plot_cdf_cell(RoA_NEG,Names_Conditions);
    xlabel('RoA'); set(gcf,'Name','NEGATIVE')
    plot_cdf_cell(ISIs_NEG,Names_Conditions);
    xlabel('ISI'); set(gcf,'Name','NEGATIVE')
    plot_cdf_cell(TranLengths_NEG,Names_Conditions);
    xlabel('ED'); set(gcf,'Name','NEGATIVE')
end