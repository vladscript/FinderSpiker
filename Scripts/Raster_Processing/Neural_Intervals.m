% Ensembles Transitions: Hebbian Pathway
% Input
%   CAGsmooth
%   frame_ensembles
%   signif_frames
% Output
%   Hebbian Sequence        Hebbian Sequence
%   Ensemble Onset          Hebbian Onsets [samples]
%   Ensemble Interval        Hebbian Intervals [samples]

function [EnsembleIntervals,EnsembleInstancesTimes,EnsembleInstances]=Neural_Intervals(ensemble_index,time_ensemble,ensemble_inter,CAGsmooth)
%% MAGIC
% CONSECUTIVE NEURAL ENSEMBLE INSTANCES ***********************************
Fsample=[find(diff(ensemble_index)~=0);numel(time_ensemble)];
Start=1;
IntervalEnsemble=[];
HebbEnsemble=[];
TimeEnsemble=[];
for f=1:numel(Fsample)
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
        % SINGLE INTERVAL
        if isempty(peaksFrames)
            IntervalEnsemble=[IntervalEnsemble;Aframe,Bframe];
            HebbEnsemble=[HebbEnsemble;ensemble_index(Start)];
            TimeEnsemble=[TimeEnsemble;round(median(time_ensemble(Start:Fsample(f))))];
        else
            IntervalEnsemble=[IntervalEnsemble;Aframe,Bframe];
            HebbEnsemble=[HebbEnsemble;ensemble_index(Start)];
            TimeEnsemble=[TimeEnsemble;round(peaksFrames(1))];
        end
    else
        % SINGLE INTERVAL during PEAK
        Npeaks=numel(peaksFrames);
        if Npeaks==1 
            IntervalEnsemble=[IntervalEnsemble;Aframe,Bframe];
            HebbEnsemble=[HebbEnsemble;ensemble_index(Start)];
            TimeEnsemble=[TimeEnsemble;round(peaksFrames(1))];
        % MULTIPLE INTERVALS
        else
            limitA=Aframe;
            for p=1:Npeaks
                if p<Npeaks
                    % Check Valley in Between
                    limitB=intersect(valleysFrames,peaksFrames(p):peaksFrames(p+1));
                else
                    limitB=Bframe;
                end
                IntervalEnsemble=[IntervalEnsemble;limitA,limitB];
                HebbEnsemble=[HebbEnsemble;ensemble_index(Start)];
                TimeEnsemble=[TimeEnsemble;round((peaksFrames(p)))];
                limitA=limitB+1;
            end
        end
    end

    Start=Fsample(f)+1;
end
%% 3th WORK-AROUND

%% UPDATE OUTPUT
% EnsembleInstancesTimes=TimesEnsembles;
% EnsembleInstances=InstancesEnsembles;
% EnsembleIntervals=get_ensemble_intervals(TimesEnsembles,InstancesEnsembles,signif_frames,labels_frames);
% EnsembleInstancesTimes=round(median(EnsembleIntervals,2));
% EnsembleInstances=InstancesEnsembles;
% Get Only Ensemables that last at least 1/2-second
% OKensInsta=find(EnsembleIntervals(:,2)-EnsembleIntervals(:,1)>round(fs/2));
% EnsembleIntervals=IntervalEnsemble;
% EnsembleInstancesTimes=TimeEnsemble;
% EnsembleInstances=HebbEnsemble;

% INTERLEAVED NEURONAL ENSEMBLES INSTANCES*********************************
[~,Fp,~,~]=findpeaks(CAGsmooth);
[~,Fv,~,~]=findpeaks(-CAGsmooth);

for n=1:numel(Fp)
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
        EnsatPeak=unique(HebbEnsemble(indxEns));
        fprintf('>>')
        % Merge Ensembles at the SAME PEAK:
        % It repears Time Instance and Intervals as redundant:
        for m=1:numel(EnsatPeak)
            % Update Ensemble Intel
            MergeIndx=find(HebbEnsemble(indxEns)==EnsatPeak(m));
            InterSameEns=IntervalEnsemble(indxEns(MergeIndx),:);
            [~,indxMaxInter]=max(InterSameEns(:,2)-InterSameEns(:,1));
            % Update Time Ensemble
            TimeEnsemble(indxEns(MergeIndx))=TimeEnsemble(indxEns(MergeIndx(indxMaxInter)));
            % Update Intervals
            NewInterval=[InterSameEns(1,1),InterSameEns(end,2)];
            IntervalEnsemble(indxEns(MergeIndx),:)=repmat(NewInterval,numel(MergeIndx),1);
            fprintf('.')
        end
        fprintf('\n')
    else
        fprintf('\n')
    end
end
%% FINAL OUTPUT ***********************************************************
[EnsembleInstancesTimes,IndexsAR]=unique(TimeEnsemble);
EnsembleIntervals=IntervalEnsemble(IndexsAR,:);
EnsembleInstances=HebbEnsemble(IndexsAR);