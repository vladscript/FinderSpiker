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
%% Set Colors
[CM,ColorIndx]=Color_Selector(Names_Conditions);
%% Plots
DodgeInit=0.3;
Dodge=DodgeInit;
aux=1;
labelsplot=[];
for c=1:NC
    for d=c:NC
        if d==c
            ActualColor=CM(c,:);
        else
            % Mix Colors:
            ActualColor=double(imadd(uint8(255*CM(c,:)),uint8(255*CM(d,:)),'uint8')/255);
        end
        hplot{c}=raincloud_plot(SIM_MATRIX{c,d},'color',ActualColor,'box_on',1,'alphaval',1,'box_dodge', 1, 'box_dodge_amount',Dodge, 'dot_dodge_amount', Dodge, 'box_col_match',0);
        labelsplot=[labelsplot,hplot{c}{1}];
        Dodge=DodgeInit*(aux+1);
        Labels{aux}=[Names_Conditions{c},' vs ',Names_Conditions{d}];
        aux=aux+1;
        % ksdensity(SimValues,linspace(min(SimValues),max(SimValues),100));
        %auxcomb=auxcomb+1;
    end
end
legend(labelsplot,Labels,'Location','northwest');
hold(ax1,'off');
axis(ax1,'tight');
ax1.XLim=[0,1];
grid(ax1,'on');

