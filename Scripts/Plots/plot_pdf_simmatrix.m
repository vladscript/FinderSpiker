% Function to plot pdf's of simmilarities
% between or among ensembles of different condition
function plot_pdf_simmatrix(SIM_MATRIX,Names_Conditions)
%% Setup
NC=numel(Names_Conditions);
auxcomb=1;
figSimmMAt=figure;
figSimmMAt.Name='PDFs for Simmilarities among ensembles';
ax1=subplot(1,1,1);
hold(ax1,'on');
for c=1:NC
    for d=c:NC
        SimValues=SIM_MATRIX{c,d};
        Labels{auxcomb}=[Names_Conditions{c},' vs ',Names_Conditions{d}];
        ksdensity(SimValues,linspace(min(SimValues),max(SimValues),100));
        auxcomb=auxcomb+1;
    end
end
hold(ax1,'off');
axis(ax1,'tight');
ax1.XLim=[0,1];
grid(ax1,'on');
legend(Labels)
