%% Figure OVAs ************************************************************
function plot_ova(x1Grid,x2Grid,posterior,CMprobs,Xpca,ONEpca,labelconditions,...
        TWOpca,Y,CM,ColorIndx,RangeProbs,Yconfident)
% setup 
Labels=unique(Y,'legacy');
numClasses=numel(Labels);

%  Plot:
fprintf('>>Making Plot: ...\n')
figure
NcolsFig=nextpow2(numClasses+1);
haxis = zeros(numClasses + 1,1); % Preallocation for graphics handles
% Ybin=getbinarylabel(Y);
for j = 1:numClasses
    fprintf('>>Surface of Model %i: ...',j)
    % MissClassifed=find(Yallbin(:,j)~=Ybin(:,j));
    Ybin= false(numel(Y),1); Ybin(Y==labelconditions(j))=true;
    MissClassifed=find(Yconfident{j}~=Ybin);
    Ax{j}=subplot(2,NcolsFig,j);
    contourf(x1Grid,x2Grid,reshape(posterior{j},size(x1Grid,1),size(x1Grid,2)));
    colormap(CMprobs);
    caxis(RangeProbs);
    hold on
    % h(1:numClasses) = gscatter(Xpca(:,ONEpca),Xpca(:,TWOpca),Y,CM(ColorIndx,:),'ooo',10);
    h = gscatter(Xpca(:,ONEpca),Xpca(:,TWOpca),Y,CM(ColorIndx,:),'ooo',8);
    % Missclassified (gray)
    gscatter(Xpca(MissClassifed,ONEpca),Xpca(MissClassifed,TWOpca),Y(MissClassifed),'kkk','xxx',12);
    for n=1:numel(ColorIndx)
        h(n).MarkerEdgeColor=[0,0,0];
        if n==j
            h(n).MarkerFaceColor=CM(ColorIndx(n),:);
        else
            % False Classes
            h(n).MarkerFaceColor=[0.8,0.8,0.8];
        end;
        h(n).LineWidth=0.5;
    end
    title(sprintf('%s vs all',char(labelconditions(j))));
    xlabel(['PC ',num2str(ONEpca)]);
    ylabel(['PC ',num2str(TWOpca)]);
    legend off
    axis tight
    hold off
    fprintf('done.\n')
end


haxis(numClasses + 1) = colorbar('Location','EastOutside',...
    'Position',[[0.8,0.1,0.05,0.4]]);
set(get(haxis(numClasses + 1),'YLabel'),'String','Posterior True Class','FontSize',9);

% legend(h(1:numClasses),'Location',[0.6,0.2,0.1,0.1]);

for j = 1:numClasses
    Axis=Ax{j};
    hlines = findobj(Axis,'Type','line');
    nlineFalseclass=0;
    nlineTrueclass=0;
    nlineMisslass=0;
    StringsCell={'True','False','Miss'};
    for n=1:numel(hlines)
        ColorActual=hlines(n).MarkerFaceColor;
        if numel(find(ColorActual==0.8))==3
            nlineFalseclass=n;
        end
        if numel(intersect( ColorActual, CM(ColorIndx(j),:)))==3
            nlineTrueclass=n;
        end
        if hlines(n).Marker=='x'
            nlineMisslass=n;
        end
    end
    legendobj=round([nlineTrueclass,nlineFalseclass,nlineMisslass]);
    legendobjsindx=find(legendobj>0);
    legend(hlines(legendobj(legendobjsindx)),StringsCell(legendobjsindx),...
        'FontSize',8,'Location','south')
end

% Legends of colors for true classes
% Legend of false classes
% 