% Funtion to re-label hebbian sequence
% Input
%   labels:     vector sequence of neural ensemble instances
%   mode: 'sequence','2-freq'
% Output
%   relabels    resorted sequence of neural ensemble instances
function relabel=relabel_ensembles(labelsEns,HebbSequence,modesort)
%% Setup
fprintf('>>[relabeling ensemble instances]: ')
relabel=zeros(size(labelsEns));
NensOK=numel(unique(labelsEns));
switch modesort
    case 'sequence'
        %% Case mode
        % RE-LABEL Ensembles
        OKsorting=unique(HebbSequence,'stable');
        
    case '2-freq'
        %% Case 2-freq
        % RE-LABEL Ensembles
        TableEns=tabulate(HebbSequence);
        [~,OKsorting]=sort(TableEns(:,3),'descend');
        DommEns=OKsorting(1);
        DommEnsMinus=OKsorting(2);
        OKsorting(1)=DommEnsMinus;
        OKsorting(2)=DommEns;
end
for n=1:NensOK
    relabel(labelsEns==OKsorting(n))=n;
end
fprintf('done.\n')