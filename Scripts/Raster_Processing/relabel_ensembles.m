% Function to re-label hebbian sequence
% Input
%   labels:         Vector sequence of neural ensemble instances
%   HebbSequence    Vector of Ensembles Sequence
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
        if numel(OKsorting)>1
            DommEnsMinus=OKsorting(2);
            OKsorting(1)=DommEnsMinus;
            OKsorting(2)=DommEns;
        end
end
for n=1:NensOK
    % Dissapeared Ensembles by labeling Sequence
    if n<=numel(OKsorting)
        relabel(labelsEns==OKsorting(n))=n;
    end
end
fprintf('done.\n')