function plot_box(RASTER_NAMES,RASTER_FEATURES,Names_Conditions,f,axis)
    NC=numel(Names_Conditions);
    % f: Feature Index
    F=[];
    COND=[];
    for c=1:NC
        F=[F;RASTER_FEATURES{c}(:,f)];
        COND=[COND;repmat(categorical(cellstr(Names_Conditions{c})),size(RASTER_FEATURES{c}(:,f)))];
    end
    boxplot(axis,F,COND,'labels',Names_Conditions)
    % 1st COndition as Reference (always)
    if NC>1
        hold(axis,'on');
        Names=RASTER_NAMES{1};
        N=numel(Names);
        for n=1:N
            for c=2:NC
                for k=1:numel(RASTER_NAMES{c})
                    Names2=cell2mat( RASTER_NAMES{c}(k,:) );
                    if strcmp(Names(n),Names2)
                        Fpre=RASTER_FEATURES{c-1}(n,f);
                        Fpost=RASTER_FEATURES{c}(k,f);
                        plot(axis,[c-1,c],[Fpre,Fpost],'k')
                    end
                end
            end
        end
    end
end