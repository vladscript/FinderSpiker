if exist('RoA_POS','var')
    aremerged=true;    % Are there already colocated cells
else
    aremerged=false;    % Are there already colocated cells
end

[CM,ColorIndx]=Color_Selector(Names_Conditions);

plot_cdf_cell(RoA_ALL,Names_Conditions,CM(ColorIndx,:));
xlabel('RoA')
plot_cdf_cell(ISIs_ALL,Names_Conditions,CM(ColorIndx,:));
xlabel('ISI')
plot_cdf_cell(TranLengths_ALL,Names_Conditions,CM(ColorIndx,:));
xlabel('ED')
plot_cdf_cell(RoT_ALL,Names_Conditions,CM(ColorIndx,:));
xlabel('Rate of Transients')
if aremerged
    % POSITIVE
    plot_cdf_cell(RoA_POS,Names_Conditions,CM(ColorIndx,:));
    xlabel('RoA'); set(gcf,'Name','POSITIVE')
    plot_cdf_cell(ISIs_POS,Names_Conditions,CM(ColorIndx,:));
    xlabel('ISI'); set(gcf,'Name','POSITIVE')
    plot_cdf_cell(TranLengths_POS,Names_Conditions,CM(ColorIndx,:));
    xlabel('ED'); set(gcf,'Name','POSITIVE')
    plot_cdf_cell(RoT_POS,Names_Conditions,CM(ColorIndx,:));
    xlabel('Rate of Transients'); set(gcf,'Name','POSITIVE')
    % NEGATIVE
    plot_cdf_cell(RoA_NEG,Names_Conditions,CM(ColorIndx,:));
    xlabel('RoA'); set(gcf,'Name','NEGATIVE')
    plot_cdf_cell(ISIs_NEG,Names_Conditions,CM(ColorIndx,:));
    xlabel('ISI'); set(gcf,'Name','NEGATIVE')
    plot_cdf_cell(TranLengths_NEG,Names_Conditions,CM(ColorIndx,:));
    xlabel('ED'); set(gcf,'Name','NEGATIVE')
    plot_cdf_cell(RoT_NEG,Names_Conditions,CM(ColorIndx,:));
    xlabel('Rate of Transients'); set(gcf,'Name','NEGATIVE')
end