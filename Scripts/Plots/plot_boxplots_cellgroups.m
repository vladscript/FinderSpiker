% Function to plot ACtivity of Cells in 3 stacks of boxplots:
% Facilitated
% Depressed
% Unchanged
function plot_boxplots_cellgroups(GlobalRoA)

figure; 

[bda_final,bdam_final]=pairedraincloud('paired',[0,1,0;0,1,0],0,0,...
    GlobalRoA.fac(:,1),GlobalRoA.fac(:,2));


[bda_final,bdam_final]=pairedraincloud('paired',[0,0,1;0,0,1],bda_final,bdam_final,...
    GlobalRoA.unc(:,1),GlobalRoA.unc(:,2));

[~,~]=pairedraincloud('paired',[1,0,0;1,0,0],bda_final,bdam_final,...
    GlobalRoA.dep(:,1),GlobalRoA.dep(:,2));

% figure
% raincloud_plot(GlobalRoA.fac(:,1),'color',CM(1,:),'box_on',1,...
%     'alphaval',5,'box_dodge',1,...
%     'box_dodge_amount',0.4, 'dot_dodge_amount', 0.4,...
%     'box_col_match',0,'box_dodge_amount',.25, 'dot_dodge_amount', 0.25,...
%     'line_width',3,'lwr_bnd',2);
% 
% raincloud_plot(GlobalRoA.fac(:,2),'color',CM(1,:),'box_on',1,...
%     'alphaval',5,'box_dodge',1,...
%     'box_dodge_amount',0.8, 'dot_dodge_amount', 0.8,...
%     'box_col_match',0,'box_dodge_amount',.5, 'dot_dodge_amount', 0.5,...
%     'line_width',3,'lwr_bnd',2);
% drawnow;
% Bar1=gca;
% h1=findobj(Bar1,'Type','Scatter');
% Xdata1=h1(1).XData;
% Ydata1=h1(1).YData;
% Xdata2=h1(2).XData;
% Ydata2=h1(2).YData;
% for i=1:numel(GlobalRoA.fac(:,2))
%     plot([Xdata1(i),Xdata2(i)],[Ydata1(i),Ydata2(i)],...
%         'LineStyle','--','Color',CM(1,:),'LineWidth',0.1)
% end