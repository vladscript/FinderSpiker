% function to plot Cummulative Distribution Function
% of the data at input cell 
function plot_cdf_cell(RoA_All,Names_Conditions,ColMap)
NC=size(ColMap,1);
CDFfig=figure;
CDFfig.Name='Cummulative Function Distribution';
ax1=subplot(1,1,1);
hold(ax1,'on');
for c=1:NC
    % if ~isempty(RoA_All{c})
    if c<=numel(RoA_All)
        if min(RoA_All{c})~=max(RoA_All{c})
            [f,xi]=ksdensity(RoA_All{c},linspace(min(RoA_All{c}),max(RoA_All{c}),100),...
                'function','cdf');
%             [f,xi]=ksdensity(RoA_All{c},'function','cdf','support',[-0.001,1]);
            plot(xi,f,'LineWidth',2,'Color',ColMap(c,:))
        else
            plot(0,0,'*');
            disp('>>NO DATA')
        end
    else
        plot(0,0,'*');
        disp('>>NO DATA')
    end
end
hold(ax1,'off');
axis(ax1,'tight');
grid(ax1,'on');
legend(Names_Conditions(1:NC))