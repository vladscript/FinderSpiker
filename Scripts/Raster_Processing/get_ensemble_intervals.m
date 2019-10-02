% Input
%   signif_frames:  all active frames
%   EnsembleTimes:  mode of ensemble instance
%   Ensembles_Labels:   ensemble index
% Output
% EnsembleInterva: starting and final time for each ensmeble instance
function EnsembleInterva=get_ensemble_intervals(EnsembleTimes,HebbSequence,signif_frames,Ensembles_Labels)
%% Setup
Nensins=numel(EnsembleTimes);
Bp=1;
EnsembleInterva=zeros(Nensins,2);
for n=1:Nensins
    if n<Nensins
        preEns=EnsembleTimes(n);
        posEns=EnsembleTimes(n+1);
        AntEns=HebbSequence(n);
        PosEns=HebbSequence(n+1);
    else
        preEns=EnsembleTimes(n);
        posEns=signif_frames(end);
        AntEns=HebbSequence(n);
        PosEns=HebbSequence(n);
    end
    preIndex=find(signif_frames>=preEns);
    posIndex=find(signif_frames<=posEns);
    preFrames=signif_frames(preIndex);
    posFrames=signif_frames(posIndex);
    
    % Check Sequence
    [InBWTN,Axx]=intersect(preFrames,posFrames);
    
    
    if AntEns==PosEns
        if n<Nensins
        % Finde where interval is maximum
        [~,nmax]=max(diff(InBWTN));
        else
            nmax=numel(InBWTN);
        end
    else
        % FInd where it is differnt ensemble
        nmax=find(diff(Ensembles_Labels(preIndex(Axx)))~=0,1);
        if isempty(nmax)
            nmax=numel(InBWTN);
        end
    end
    
    nA=find(posFrames>max(sort(EnsembleInterva(:))),1);
    Ap=posFrames(nA);

    Bp=InBWTN(nmax);
    if Ensembles_Labels(signif_frames==Bp)~=AntEns
        LabeInBwtn=[];
        auxj=1;
        for j=Ap:InBWTN(end)
            if ~isempty(Ensembles_Labels(signif_frames==j));
                LabeInBwtn(1,auxj)=Ensembles_Labels(signif_frames==j);
                FrameLab(1,auxj)=signif_frames(signif_frames==j);
                auxj=auxj+1;
            end
            
        end
        oklabelindx=find(LabeInBwtn==AntEns);
        Bp=FrameLab(oklabelindx(end));
    end
    
    EnsembleInterva(n,1)=Ap;
    EnsembleInterva(n,2)=Bp;
end