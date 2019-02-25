% Plot groups in raster
% Plot the different groups in raster.
%
% PlotGroupsR(X, PeaksIdx, PeaksGroups, numFig)
%
% Inputs
%
%
% Output
% 
%
% ..:: by Jesús E. Pérez-Ortega ::.. Jun-2012 Agu-2012

function PlotGroupsR_JP(X, PeaksIdx, PeaksGroups, g2plot, numFig)

ptsymb = {'.r','.g','.b','.c','.m','.y','.k'};
[F C]=size(X);
Co=sum(X,2);
peaks=max(PeaksIdx);

figure(numFig);
subplot(5,5,[1:4 6:9 11:14 16:19])
hold on
for i=1:peaks
    if find(g2plot==PeaksGroups(i))
        peak=find(PeaksIdx==i);
        for j=1:length(peak)
            idxActive=find(X(peak(j),:));
            plot(peak(j),idxActive,ptsymb{PeaksGroups(i)},'MarkerSize',10)
        end
    end
end
hold off

subplot(5,5,[21:24])
hold on
for i=1:peaks
    if find(g2plot==PeaksGroups(i))
        peak=find(PeaksIdx==i);
        plot(peak,Co(peak),ptsymb{PeaksGroups(i)})
    end
end