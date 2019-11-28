% Ensembles Transitions: Hebbian Pathway
% Input
%   CAGsmooth
%   frame_ensembles
%   signif_frames
% Output
%   Hebbian Sequence        Hebbian Sequence
%   Ensemble Onset          Hebbian Onsets [samples]
%   Ensemble Interval       Hebbian Intervals [samples]

function [EnsembleIntervals,EnsembleInstancesTimes,EnsembleInstances]=Neural_Intervals(ensemble_index,time_ensemble,ensemble_inter,CAGsmooth)
%% CONSECUTIVE NEURAL ENSEMBLE INSTANCES **********************************
fprintf('\n>>Checking Consecutive Neuronal Ensemble Instances.\n')
Fsample=[find(diff(ensemble_index)~=0);numel(time_ensemble)];
Start=1;
IntervalEnsemble=[];
HebbEnsemble=[];
TimeEnsemble=[];
if ~isempty(ensemble_index)
    for f=1:numel(Fsample)
        fprintf('>>Neuronal Ensemble %i  at Burst %i/%i:',ensemble_index(Start),f,numel(Fsample))
        EnsembleIntervalsWhole=ensemble_inter(Start:Fsample(f),:);
        Aframe=EnsembleIntervalsWhole(1,1);
        Bframe=EnsembleIntervalsWhole(end,2);
        if Bframe-Aframe>3
            % Peaks
            [~,peaksFrames]=findpeaks( CAGsmooth(Aframe:Bframe) );
            % Valleys
            [~,valleysFrames]=findpeaks( -CAGsmooth(Aframe:Bframe) );
            peaksFrames=peaksFrames+Aframe-1;
            valleysFrames=valleysFrames+Aframe-1;
        else
            peaksFrames=[];
            valleysFrames=[];
        end
        % Intervals:
        if or(isempty(peaksFrames),isempty(valleysFrames))
            % SINGLE INTERVAL: VALLEY and NO PEAK
            if isempty(peaksFrames) && ~isempty(valleysFrames)
                limitA=Aframe;
                valleysFrames(end+1)=Bframe;
                Nvalleys=numel(valleysFrames);
                for v=1:Nvalleys
                    LimitsB=EnsembleIntervalsWhole(EnsembleIntervalsWhole(:,2)<=valleysFrames(v),2);
                    if isempty(LimitsB)
                        LimitsB=Bframe;
                    end
                    limitB=LimitsB(end);
                    IntervalEnsemble=[IntervalEnsemble;limitA,limitB];
                    HebbEnsemble=[HebbEnsemble;ensemble_index(Start)];
                    [~,MaxCAG]=max(CAGsmooth(limitA:limitB));
                    TimeEnsemble=[TimeEnsemble; MaxCAG+limitA-1 ];
                    
                    limitsA=EnsembleIntervalsWhole(EnsembleIntervalsWhole(:,1)>limitB,1);
                    if isempty(limitsA)
                        limitsA=Aframe;
                    end
                    if v<Nvalleys
                        limitA=limitsA(1);
                    end
                    fprintf('*')
                end
            % SINGLE INTERVAL: NO VALLEY and NO PEAK
            else
                IntervalEnsemble=[IntervalEnsemble;Aframe,Bframe];
                HebbEnsemble=[HebbEnsemble;ensemble_index(Start)];
                [~,MaxCAG]=max(CAGsmooth(Aframe:Bframe));
                TimeEnsemble=[TimeEnsemble; MaxCAG+Aframe-1 ];
                fprintf('*')
            end
        else
            % SINGLE INTERVAL during PEAK and VALLEY(S)
            Npeaks=numel(peaksFrames);
            if Npeaks==1
                limitA=Aframe;
                valleysFrames(end+1)=Bframe;
                Nvalleys=numel(valleysFrames);
                for v=1:Nvalleys
                    LimitsB=EnsembleIntervalsWhole(EnsembleIntervalsWhole(:,2)<=valleysFrames(v),2);
                    if isempty(LimitsB)
                        LimitsB=Bframe;
                    end
                    limitB=LimitsB(end);
                    IntervalEnsemble=[IntervalEnsemble;limitA,limitB];
                    HebbEnsemble=[HebbEnsemble;ensemble_index(Start)];
                    if ismember(peaksFrames(1),limitA:limitB)
                        TimeEnsemble=[TimeEnsemble;round(peaksFrames(1))];
                    else
                        [~,MaxCAG]=max(CAGsmooth(limitA:limitB));
                        TimeEnsemble=[TimeEnsemble; MaxCAG+limitA-1 ];
                    end
                    
                    limitsA=EnsembleIntervalsWhole(EnsembleIntervalsWhole(:,1)>limitB,1);
                    if isempty(limitsA)
                        limitsA=Aframe;
                    end
                    if v<Nvalleys
                        limitA=limitsA(1);
                    end
                    fprintf('*')
                end
            % MULTIPLE INTERVALS
            else
                % Check Peaks, Valleys and Intervals:
                if and(sum(EnsembleIntervalsWhole(:)>valleysFrames(end))>0,...
                    valleysFrames(end)>peaksFrames(end))
                    peaksFrames(end+1)=Bframe;
                end
                Npeaks=numel(peaksFrames);
                limitA=Aframe;
                for p=1:Npeaks
                    if p<Npeaks
                        % Check Valley in Between Actual and Next Peak
                        valleyB=intersect(valleysFrames,peaksFrames(p):peaksFrames(p+1));
                        LimitsB=EnsembleIntervalsWhole(EnsembleIntervalsWhole(:,2)<=valleyB,2);
                        limitB=LimitsB(end);
                    else
                        valleyB=intersect(valleysFrames,peaksFrames(p):Bframe);
                        if ~isempty(valleyB)
                            LimitsB=EnsembleIntervalsWhole(EnsembleIntervalsWhole(:,2)<=valleyB,2);
                            limitB=LimitsB(end);
                        else
                            limitB=Bframe;
                        end
                    end
                    % To ignore minor peaks in between Limits
                    if limitA<=limitB
                        % OUTPUT
                        IntervalEnsemble=[IntervalEnsemble;limitA,limitB];
                        HebbEnsemble=[HebbEnsemble;ensemble_index(Start)];
                        if ismember(peaksFrames(p),limitA:limitB)
                            TimeEnsemble=[TimeEnsemble;round(peaksFrames(p))];
                        else
                            [~,MaxCAG]=max(CAGsmooth(limitA:limitB));
                            TimeEnsemble=[TimeEnsemble; MaxCAG+limitA-1 ];
                        end
                        fprintf('*')
                        limitsA=EnsembleIntervalsWhole(EnsembleIntervalsWhole(:,1)>limitB,1);
                        if p<Npeaks
                            limitA=limitsA(1);
                        end
                    end
                end
            end
        end
        Start=Fsample(f)+1;
        fprintf('\n')
    end
end
fprintf('>>Checking Consecutive Neuronal Ensemble Instances Done.\n')
%% INTERLEAVED NEURONAL ENSEMBLES INSTANCES********************************
fprintf('>>Checking Interleaved Neuronal Ensemble Instances.\n')
% From the smoothed CAG check if Ensemble Instances belong to CAG Peaks
[~,Fp,~,~]=findpeaks(CAGsmooth);
[~,Fv,~,~]=findpeaks(-CAGsmooth);

for n=1:numel(Fp)
    fprintf('>>Checking Peak %i/%i\n',n,numel(Fp));
    % Check Previous and Next Valley of curretn Peak (if any)
    IndxBefore=find(Fv<Fp(n));
    IndxAfter=find(Fv>Fp(n));
    if isempty(IndxBefore)
        A=1;
    else
        A=Fv(IndxBefore(end));
    end
    if isempty(IndxAfter)
        B=numel(CAGsmooth);
    else
        B=Fv(IndxAfter(1));
    end
    % PeakInterval(n,:)=[A,B];
    [~,indxEns]=intersect( TimeEnsemble, A:B);
    if ~isempty(indxEns)
        fprintf('Neuronal Ensemble: ');
        % Enemble
        EnsatPeak=unique(HebbEnsemble(indxEns));
        % Merge Ensembles at the SAME PEAK:
        % It repears Time Instance and Intervals as redundant:
        for m=1:numel(EnsatPeak)
            % Update Ensemble Intel
            fprintf('%i,',EnsatPeak(m))
            MergeIndx=find(HebbEnsemble(indxEns)==EnsatPeak(m));
            InterSameEns=IntervalEnsemble(indxEns(MergeIndx),:);
            [~,indxMaxInter]=max(InterSameEns(:,2)-InterSameEns(:,1));
            % Update Time Ensemble
            TimeEnsemble(indxEns(MergeIndx))=TimeEnsemble(indxEns(MergeIndx(indxMaxInter)));
            % Update Intervals
            NewInterval=[InterSameEns(1,1),InterSameEns(end,2)];
            IntervalEnsemble(indxEns(MergeIndx),:)=repmat(NewInterval,numel(MergeIndx),1);
        end
    else
        fprintf('\n')
    end
    fprintf('\n')
end
%% FINAL OUTPUT ***********************************************************
[EnsembleInstancesTimes,IndexsAR]=unique(TimeEnsemble);
EnsembleIntervals=IntervalEnsemble(IndexsAR,:);
EnsembleInstances=HebbEnsemble(IndexsAR);