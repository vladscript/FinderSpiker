%% Script to plot Hebbian Stairs: *************************************
% Sequence of neuronal ensembles
%   Ouput Figure of same axis limits
% N conditions
% Hebbian Sequence (vector)
% Color of Ensembles
%% #####################################################################
if exist('Features_Condition','var')
    EnsembleSequence=[];
    N=numel(Features_Condition.HebbianSeq);
    Naux=1;
    LimitX=0;
    for n=1:N
        figure('NumberTitle','off','Name',['Neuronal Ensembles Sequence ',Names_Conditions{n}],...
            'Position',[730 470-40*(n-1) 560 151]); hold on;
        % LINE ------------------------------------------------------------
        plot(Naux:Naux+numel(Features_Condition.HebbianSeq{n})-1,...
            Features_Condition.HebbianSeq{n},...
        'LineWidth',2,'Color','k');
        Nens=numel(unique(EnsembleSequence));
        % CIRCLES ooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
        aux_ne=1;
        for n_e=Naux:Naux+numel(Features_Condition.HebbianSeq{n})-1
            plot(n_e,Features_Condition.HebbianSeq{n}(aux_ne),'o',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',ColorState(Features_Condition.HebbianSeq{n}(aux_ne),:),...
            'MarkerSize',15);
            aux_ne=aux_ne+1;
        end
        Naux=numel(Features_Condition.HebbianSeq{n})+1;
        EnsembleSequence=[EnsembleSequence;Features_Condition.HebbianSeq{n}];
        hold off;
        ylabel('Neuronal Ensemble')
        xlabel('Neuronal Ensemble Realization')
        CurrAxis{n}=gca;
        Ydiff(n)=diff(CurrAxis{n}.YLim);
        CurrAxis{n}.XLim(1)=LimitX;
        LimitX=CurrAxis{n}.XLim(end);
    end
    [MaxAxisY,INdexY]=max(Ydiff);
    for n=1:N
        CurrAxis{n}.YLim(2)=CurrAxis{n}.YLim(2)+MaxAxisY-diff(CurrAxis{n}.YLim);
    end
    fprintf('####################################################\n')
    fprintf('#                                                  #\n')
    fprintf('#  Check Y- and X- axis match for more experiments #\n')
    fprintf('#                                                  #\n')
    fprintf('####################################################\n\n')
else
    fprintf('Missing Clustering Analysis, please execute:\n')
    fprintf('>>R_CONDITIONi_Analysis=get_bays_ensembles(R_CONDTIONi);\n')
    fprintf('For every Raster Condition\n')
end
%% END