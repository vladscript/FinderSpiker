function FinalTable=metrotable(METRICS,labelconditions)
FinalTable=table;
for c=1:numel(METRICS)
    % Multiclass Binary Complete * * * * * * * * * * * * * * * * * * * * *
    MCmetrics=METRICS{c};
    MetricsNames=fieldnames(MCmetrics);
    RowName{c}=char(labelconditions(c));
    RowString=cell(1,numel(MetricsNames));
    for n=1:numel(MetricsNames)-1
        % Columna Name
        metric=getfield(MCmetrics, MetricsNames{n});
        RowString{n}=sprintf('%3.2f(%3.2f)',mean(metric),std(metric));
    end
    RowTable=RowString;
    FinalTable=[FinalTable;RowTable];
end
FinalTable.Properties.RowNames=RowName;
FinalTable.Properties.VariableNames=MetricsNames;