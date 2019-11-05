% Ensembles Transitions: Hebbian Pathway
% Input
%   fs:             sampling frequency (4 plot only)
%   labels_frames:  Ensemble Lalbel
%   signif_frames:  Frames with Significative CoACtivity
%   CAG:            Coactivitygraphy
%   ColorState:     Ensemble Color      [Nensembles,3]
%   If              ifplot=1 plot; else do nothing
%   Limist of COnditions      [frames] NC>1
% Output
%   EnsembleInstances:      Hebbian Sequence
%   EnsembleInstancesTimes: Hebbian Onsets [samples]
%   EnsembleIntervals:      Hebbian Duration [samples]
function [EnsembleInstances,EnsembleInstancesTimes,EnsembleIntervals]=Ensembles_Transitions(fs,labels_frames,signif_frames,CAG,ColorState,ifplot,varargin)
%% Setup
% Ignore 0-ensembled frames
Load_Default_Clustering;
EnsembleInstances=[];
EnsembleIntervals=[];
% EnsemblesAllIntervals=[];
labfram=labels_frames(labels_frames>0);
framsig=signif_frames(labels_frames>0);
frames_ensembles=[];
aux=1;
time_ensemble=[];
ensemble_index=[];
HebbianSequence=[];
TimesHebb=[];
labfram=[labfram;-1000];
framsig=[makerowvector(framsig),1000];

%% Read Continuous Ensembled-Frames
for i=1:numel(framsig)-1;
    % Find Consecutive Frames
    if framsig(i+1)==framsig(i)+1 &&...
        labfram(i)==labfram(i+1)
        frames_ensembles=[frames_ensembles,framsig(i)];
    else
        if isempty(frames_ensembles)
            frames_ensembles=framsig(i);
        end
        %  --------------------------------------------
        time_ensemble(aux,1)=round(mean(frames_ensembles)); % $$$
        ensemble_index(aux,1)=labfram(i);                   % $$$
        aux=aux+1;
        frames_ensembles=[];
    end 
end

% Smooth Hebbian Sequence:
% Th=min(CAG(time_ensemble));
%% IF CAG intel::
if isempty(CAG)
    EnsembleInstancesTimes=time_ensemble;
    EnsembleInstances=ensemble_index;
    % EnsembleIntervals=[time_ensemble-1,time_ensemble+1];
    % Every ball is an ensmble of continous frames
else
    % Smoothed Version of CAG: Ensmble above CAG Threshold
    % Read Frame Limits and CAG Thresholds:
    NC=numel(varargin{1}); % N conditions
    % Get zero-CAGs
    %ZeroCAGAll=find(CAG==0);
    % Ignore one-second-frame Zeros
    XeroDrive=zeros(size(CAG));
    XeroDrive(CAG==0)=1;
    dXero=diff(XeroDrive);
    StartsIntervasl=find(dXero>0);
    EndsIntervasl=find(dXero<0);
    ZeroCAG=[];
    if and(~isempty(StartsIntervasl),~isempty(EndsIntervasl))
        % If it started at Initial Transient
        if StartsIntervasl(1)>EndsIntervasl(1)
            %EndsIntervasl=EndsIntervasl(2:end);
            StartsIntervasl=[0,makerowvector(StartsIntervasl)];
        end
        % If it ended at with out End of Transient
        if StartsIntervasl(end)>EndsIntervasl(end)
            % StartsIntervasl=StartsIntervasl(1:end-1);
            EndsIntervasl=[makerowvector(EndsIntervasl),numel(CAG)];
        end
        
        if numel(StartsIntervasl)==numel(EndsIntervasl)
            for f=1:numel(StartsIntervasl)
                if EndsIntervasl(f)-StartsIntervasl(f)>fs*CAGDeepValleys
                    ZeroCAG=[ZeroCAG,StartsIntervasl(f)+1:EndsIntervasl(f)];
                end
            end
        else
            disp('WTF?');
        end
    end
    if NC==1
        fprintf('>>One Condition Detected')
        % CAG Threshold 
        CAGth=min(CAG(signif_frames));
        % Smoothed CAG
        Nsamp2smoot=round(numel(CAG)/40);
        if Nsamp2smoot<2 
            Nsamp2smoot=2; 
        end
        CAGsmooth=smooth(CAG,Nsamp2smoot,'loess');
        % Deep Valleys
        [~,ValleysF]=findpeaks(-CAGsmooth);
        DeepValleys=ValleysF(CAGsmooth(ValleysF)<CAGth);
    else
        fprintf('>>Several Conditions Detected')
        % Initiate
        LIMS=varargin{1};
        Af=1;
        CAGsmooth=zeros(size(CAG));
        ValleysF=[];
        DeepValleys=[];
        for c=1:NC
            Bf=LIMS{c}+Af-1;
            FramConds=intersect(signif_frames,Af:Bf);
            % CAG Thresholds
            if ~isempty(FramConds)
                CAGth(c)=min(CAG(FramConds));
            else
                CAGth(c)=0;
            end
            % Smoothed CAGs
            Nsamp2smoot=round(numel(CAG(Af:Bf))/40);
            if Nsamp2smoot<2 
                Nsamp2smoot=2; 
            end
            CAGsmooth(Af:Bf)=smooth(CAG(Af:Bf),Nsamp2smoot,'loess');
            CAGc=CAGsmooth(Af:Bf); % Condition vector
            % Deep Valleys
            [~,ValleysFc]=findpeaks(-CAGc);
            ValleysF=[ValleysF;ValleysFc+Af-1];
            if CAGth(c)>0
                DeepValleysC=ValleysFc(CAGc(ValleysFc)<CAGth(c));
            else
                DeepValleysC=[];
            end
            DeepValleys=[DeepValleys;DeepValleysC+Af-1];
            Af=Bf+1;
        end
    end
    ValleysF=[1;ValleysF;numel(CAG)];
    %% 1st WORK-AROUND
    
    Fsample=[find(diff(ensemble_index)~=0);numel(time_ensemble)];
    Start=1;
    for f=1:numel(Fsample)
        EnsembleInterval=time_ensemble(Start:Fsample(f));
        % Check Vallleys in the Interval to determine Smoothed Sequence
        PrevAlleys=find(ValleysF<=EnsembleInterval(1));
        PostAlleys=find(ValleysF>=EnsembleInterval(end));
        if ~isempty( PrevAlleys )
            prevalley=PrevAlleys(end);
        else
            prevalley=1;
        end
        if ~isempty(PostAlleys)
            posvalley=PostAlleys(1);
        else
            posvalley=numel(PostAlleys);
        end

        NCircles=posvalley-prevalley;
        for n=1:NCircles
            if NCircles==numel(EnsembleInterval)
                ActualHebb=EnsembleInterval(n);
                ActualEnsemble=ensemble_index(Start);
            else
                TimeA=find(EnsembleInterval>ValleysF(prevalley),1);
                TimeB=find(EnsembleInterval<ValleysF(prevalley+1));
                TimeB=TimeB(end);
                if TimeB>TimeA
                    ActualHebb=round(median(EnsembleInterval(TimeA:TimeB)));
                    ActualEnsemble=ensemble_index(Start);
                elseif TimeB==TimeA
                    ActualHebb=round(median(EnsembleInterval(TimeB)));
                    ActualEnsemble=ensemble_index(Start);
                else
                    ActualHebb=[];
                    ActualEnsemble=[];
                end
            end
            % Repeated Ensemble's Instances
            if ismember(ActualHebb,TimesHebb)
                ActualHebb=[];
                ActualEnsemble=[];
            end
            % [~,MaxLocal]=max(CAGsmooth(ValleysF(prevalley):ValleysF(prevalley+1)));
            % ActualHebb=ValleysF(prevalley)+MaxLocal-1;
            %if ~ismember(ActualHebb,TimesHebb)
            TimesHebb=[TimesHebb,ActualHebb];
            % else
            % disp('')
            %    ActualHebb=round(mean(ActualHebb:ValleysF(prevalley+1)));
            %    TimesHebb=[TimesHebb,ActualHebb];
            %end 
            HebbianSequence=[HebbianSequence,ActualEnsemble];
            prevalley=prevalley+1;
        end
        Start=Fsample(f)+1;
    end
    %% 2nd WORK-AROUND
    % Ignore Faste Changing Sequences
   
    TimHebbian=[];
    EnsHebbian=[];
    EnsList=unique(HebbianSequence);
    for n=1:numel(unique(HebbianSequence))
        fprintf('\n>> [Accounting Ensemble Instances: %i]: ',EnsList(n))
        TimesHebbE=TimesHebb(HebbianSequence==EnsList(n));
        %TimesHebbE()=CAGsmooth(TimesHebb(HebbianSequence==n));
        ActualTime=TimesHebbE(1);
        Nhebbs=1;
        Nseq=1;
        EnsInt={};
        LastValley=find(DeepValleys>ActualTime,1);
        DeepValleyLim=DeepValleys(LastValley);
        DeepValleyPre=DeepValleyLim;
        BufferEnsemble=[];
        while Nhebbs<=numel(TimesHebbE)
            LastValley=find(DeepValleys>ActualTime,1);
            DeepValleyLim=DeepValleys(LastValley);
            % Re-define Sequence: if it's the same limit
            if DeepValleyPre==DeepValleyLim
                % Accumulate Ensembles
                BufferEnsemble=[BufferEnsemble;ActualTime];
                EnsInt{Nseq}=BufferEnsemble;
            else
                Nseq=Nseq+1;
                % Re-save Sequence
                BufferEnsemble=[];
                BufferEnsemble=[BufferEnsemble;ActualTime];
                EnsInt{Nseq}=BufferEnsemble;
                % Update Limit
                DeepValleyPre=DeepValleyLim;
                fprintf('-')
            end

            Nhebbs=Nhebbs+1;
            if Nhebbs<=numel(TimesHebbE)
                ActualTime=TimesHebbE(Nhebbs);
            end
            fprintf('*')
        end
        % OUTPUT: Get Ensemble Time & Ensemble
        for s=1:Nseq
            % TimHebbian=[TimHebbian;round(median(EnsInt{s}))];
            % Prevent from repeated Ensemble Onsets
            if ~ismember(round(median(EnsInt{s})),TimHebbian)
                TimHebbian=[TimHebbian;round(median(EnsInt{s}))];
            else
                if numel(EnsInt{s}(end))>1
                    TimHebbian=[TimHebbian;round(EnsInt{s}(end))];
                else
                    TimHebbian=[TimHebbian;round(EnsInt{s}(end)+1)];
                end
            end
            
            EnsHebbian=[EnsHebbian;EnsList(n)];
        end
        fprintf('\n')
    end
    [EnsembleInstancesTimes,SortingIndx]=sort(TimHebbian);
    EnsembleInstances=EnsHebbian(SortingIndx);
    %% 3th WORK-AROUND
    % Only One Ensemble Instance by Interval
    BufferTimes=[];
    TimesEnsembles=[];
    InstancesEnsembles=[];
    % Add last unpossible ensemble & Time:
    EnsembleInstances=[EnsembleInstances;-666];
    EnsembleInstancesTimes=[EnsembleInstancesTimes;0];
    if numel(EnsembleInstances)>1
        preEns=EnsembleInstances(1);
        for n=1:numel(EnsembleInstances)
            if EnsembleInstances(n)==preEns
                % The same Ensemble
                BufferTimes=[BufferTimes;EnsembleInstancesTimes(n)];
            else
                % Ca$h cut
                if numel(BufferTimes)>1
                    % Check if CAG touches Zeros
                    prenn=1;
                    Sequ=1;
                    NNindx={};
                    nnIndexes=[];
                    for nn=1:numel(BufferTimes)
                        EnsInterval=BufferTimes(prenn):BufferTimes(nn);
                        ZeroEns=intersect(EnsInterval, ZeroCAG );
                        if isempty(ZeroEns)
                            % Accumulate
                            nnIndexes=[nnIndexes,prenn,nn];
                        else
                            % Separate
                            NNindx{Sequ}=unique(nnIndexes);
                            Sequ=Sequ+1;
                            nnIndexes=[];
                            nnIndexes=[nnIndexes,nn];
                        end
                        prenn=nn;
                    end
                    if isempty(NNindx)
                        TimesEnsembles=[TimesEnsembles;round(median(BufferTimes))];
                        InstancesEnsembles=[InstancesEnsembles;preEns];
                        % EnsemblesAllIntervals=[EnsemblesAllIntervals;BufferTimes(end)-BufferTimes(1)+1];
                    else
                        
                        if NNindx{end}(end)<numel(BufferTimes)
                            NNindx{end+1}=numel(BufferTimes);
                        end
                        for nn=1:numel(NNindx)
                            TimesEnsembles=[TimesEnsembles;round(median(BufferTimes(NNindx{nn})))];
                            InstancesEnsembles=[InstancesEnsembles;preEns];
                            % EnsemblesAllIntervals=[EnsemblesAllIntervals;BufferTimes(end)-BufferTimes(1)+1];
                        end
                    end
                else
                    TimesEnsembles=[TimesEnsembles;round(median(BufferTimes))];
                    InstancesEnsembles=[InstancesEnsembles;preEns];
                    % EnsemblesAllIntervals=[EnsemblesAllIntervals;BufferTimes(end)-BufferTimes(1)+1];
                end
                % Different Ensemble
                preEns=EnsembleInstances(n);
                BufferTimes=[];
                BufferTimes=[BufferTimes;EnsembleInstancesTimes(n)];
            end
        end
    else
        fprintf('>>Low Neural Ensemble Activity')
    end
    %% UPDATE OUTPUT
    % EnsembleInstancesTimes=TimesEnsembles;
    % EnsembleInstances=InstancesEnsembles;
    EnsembleIntervals=get_ensemble_intervals(TimesEnsembles,InstancesEnsembles,signif_frames,labels_frames);
    EnsembleInstancesTimes=round(median(EnsembleIntervals,2));
    EnsembleInstances=InstancesEnsembles;
    % Get Only Ensemables that last at least 1/2-second
    OKensInsta=find(EnsembleIntervals(:,2)-EnsembleIntervals(:,1)>round(fs/2));
    EnsembleIntervals=EnsembleIntervals(OKensInsta,:);
    EnsembleInstancesTimes=EnsembleInstancesTimes(OKensInsta);
    EnsembleInstances=EnsembleInstances(OKensInsta);
end

%% FIGURE *****************************************************************
% Requires previous figure with Raster Plot
if ifplot
    %% Plot in Minutes
    % Always plot after a Raster Analysis to get the axis right
    Axis_details=gca;
    
    
    % PLTO INTERVALS AND COLORS
    
    xtime=Axis_details.Children(end).XData; % [MINUtES]
    yCAG=Axis_details.Children(end).YData;
    
    
    fig_Ensembles_Transitions=figure('Position', [415 172 560 122],...
        'Name','Hebbian Pathways');
    plot(EnsembleInstancesTimes/fs/60,EnsembleInstances,'k','LineWidth',2); hold on
    Ntran=length(EnsembleInstances);
    for i=1:Ntran
        if ~isempty(CAG)
            % PLOT @ CAG
            if diff([EnsembleIntervals(i,1),EnsembleIntervals(i,2)])~=0
                plot(Axis_details,xtime(EnsembleIntervals(i,1):EnsembleIntervals(i,2)),...
                yCAG(EnsembleIntervals(i,1):EnsembleIntervals(i,2)),...
                'Color',ColorState(EnsembleInstances(i),:) );
            else
                plot(Axis_details,xtime(EnsembleIntervals(i,1)-1:EnsembleIntervals(i,2)+1),...
                yCAG(EnsembleIntervals(i,1)-1:EnsembleIntervals(i,2)+1),...
                'Color',ColorState(EnsembleInstances(i),:) );
            end
        end
        % PLOT @ HEBBIAN SEQUENCE
        plot(EnsembleInstancesTimes(i)/fs/60,EnsembleInstances(i),'o',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',ColorState(EnsembleInstances(i),:),...
            'MarkerSize',10); 
    end
    hold off;
    axis([Axis_details.XLim,0.5,max(EnsembleInstances)+0.5])
    EnsmbleAxis=gca;
    linkaxes([Axis_details,EnsmbleAxis],'x')
end