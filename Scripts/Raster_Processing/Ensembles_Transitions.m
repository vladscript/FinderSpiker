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
function [EnsembleInstances,EnsembleInstancesTimes,EnsembleIntervals]=Ensembles_Transitions(fs,labels_frames,signif_frames,CAG,ColorState,ifplot,R,varargin)
%% Setup
% Ignore 0-ensembled frames
Load_Default_Clustering; % See Details.
SmoothWin=SmoothWinSec*fs;
EnsembleIntervals=[];
labfram=labels_frames(labels_frames>0);
framsig=signif_frames(labels_frames>0);
frames_ensembles=[];
aux=1;
time_ensemble=[];
ensemble_index=[];
labfram=[labfram;-1000];
framsig=[makerowvector(framsig),1000];
% Neurons Ratio Proportion
EnsemblesIndex=unique(labels_frames);
if ~isempty(R)
    for ne=1:numel(EnsemblesIndex)
        CellatEns(ne)= numel(find(sum(R(:,signif_frames( labels_frames==EnsemblesIndex(ne) )),2)))/size(R,1);
    end
else
    CellatEns(1:numel(EnsemblesIndex))=1;
end
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
        time_ensemble(aux,1)=round(median(frames_ensembles));   % $$$
        ensemble_index(aux,1)=labfram(i);                       % $$$
        ensemble_inter(aux,1)=frames_ensembles(1);              % $$$
        ensemble_inter(aux,2)=frames_ensembles(end);            % $$$
        aux=aux+1;
        frames_ensembles=[];
    end 
end
if isempty(time_ensemble)
    ensemble_inter=[];
end
%% IF CAG intel::
if isempty(CAG)
    EnsembleInstancesTimes=time_ensemble;
    EnsembleInstances=ensemble_index;    
    % Every ball is an ensmble of continous frames
else
    % Smoothed Version of CAG: Ensmble above CAG Threshold
    % Read Frame Limits and CAG Thresholds:
    NC=numel(varargin{1}); % N conditions
    
    if NC==1
        fprintf('>>One Condition Detected')
        % CAG Threshold 
        CAGth=min(CAG(signif_frames));
        % Smoothed CAG
        CAGsmooth=smooth(CAG,SmoothWin,'loess');
    else
        fprintf('>>Several Conditions Detected')
        % Initiate
        LIMS=varargin{1};
        Af=1;
        CAGsmooth=zeros(size(CAG));
%         ValleysF=[];
%         Promin=[];
%         DeepValleys=[];
        for c=1:NC
            Bf=LIMS{c}+Af-1;
            FramConds=intersect(signif_frames,Af:Bf);
            % CAG Thresholds
            if ~isempty(FramConds)
                CAGth(c)=min(CAG(FramConds));
            else
                CAGth(c)=0;
            end
%             % Smoothed CAGs
            CAGsmooth(Af:Bf)=smooth(CAG(Af:Bf),SmoothWin,'loess');
            CAGc=CAGsmooth(Af:Bf); % Condition vector
            Af=Bf+1;
        end
    end

    %% Get Temporal Features of Neural Ensembles Activations---------------
    [EnsembleIntervals,EnsembleInstancesTimes,EnsembleInstances]=Neural_Intervals(ensemble_index,time_ensemble,ensemble_inter,CAGsmooth);
    % Filtering - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if ~isempty(EnsembleInstances)
        fprintf('>>Filtering Ensembles Time Features: ')
        DurEn=EnsembleIntervals(:,2)-EnsembleIntervals(:,1);
        IgnoreInstances=[];
        for n=2:numel(EnsembleInstances)
            % If interfere consecutive & too short
            if and(EnsembleIntervals(n,1)-EnsembleIntervals(n-1,2)<0,...
                    DurEn(n)<fs)
                IgnoreInstances=[IgnoreInstances;n];
                fprintf('o')
            end
        end
        OKinstances=setdiff(1:numel(EnsembleInstances),IgnoreInstances);
        % UPDATE INTEL
        EnsembleIntervals=EnsembleIntervals(OKinstances,:);
        EnsembleInstancesTimes=EnsembleInstancesTimes(OKinstances);
        EnsembleInstances=EnsembleInstances(OKinstances);
        fprintf('ok\n')
    end
end

%% FIGURE *****************************************************************
% Requires previous figure with Raster Plot
if ifplot
    %% Plot in Minutes
    % Always plot after a Raster Analysis to get the axis right
    Axis_details=gca;
    % PLOT INTERVALS AND COLORS    
    xtime=Axis_details.Children(end).XData; % [MINUtES]
    yCAG=Axis_details.Children(end).YData;
     % Figure Object for the Hebbian Sequence
    fig_Ensembles_Transitions=figure('Position', [403 172 560 122],...
        'NumberTitle','off',...
        'Name','Hebbian Sequence & Neuronal Ensmble Intervals');
    axis; HebbSeq=gca;
    if ~isempty(CAG)
        % Smoothed CAG
        plot(Axis_details,xtime,CAGsmooth,'Color',[0.9,0.9,0.9],'LineWidth',2);
        % plot(EnsembleInstancesTimes/fs/60,EnsembleInstances,'k','LineWidth',0.5); 
        [xseq,yseq]=makelineseq(EnsembleInstancesTimes,EnsembleInstances,EnsembleIntervals);
        % Line Between Neuronal
        plot(HebbSeq,xtime(xseq),yseq,'k','LineWidth',0.5);
    end
   
    hold on
    Ntran=numel(EnsembleInstances);
    fprintf('>>Hebbian Sequence:\n')
    if ~isempty(CAG)
        for i=1:Ntran
            fprintf('Instance %i - Ensemble %i - Duration %3.1f s\n',i,...
            EnsembleInstances(i),(1+EnsembleIntervals(i,2)-EnsembleIntervals(i,1))/fs);

            % PLOT @ CAG [RASTER FIGURE]
            if diff([EnsembleIntervals(i,1),EnsembleIntervals(i,2)])~=0
                plot(Axis_details,xtime(EnsembleIntervals(i,1):EnsembleIntervals(i,2)),...
                yCAG(EnsembleIntervals(i,1):EnsembleIntervals(i,2)),...
                'Color',ColorState(EnsembleInstances(i),:),'LineWidth',1);
            else
                plot(Axis_details,xtime(EnsembleIntervals(i,1)-1:EnsembleIntervals(i,2)+1),...
                yCAG(EnsembleIntervals(i,1)-1:EnsembleIntervals(i,2)+1),...
                'Color',ColorState(EnsembleInstances(i),:),'LineWidth',1);
            end
        % Plot LINE Intervals @ HEBBIAN SEQUENCE
        [~,IndxDixk]=intersect(xseq,EnsembleIntervals(i,1):EnsembleIntervals(i,2));
        plot(xtime(EnsembleIntervals(i,1):EnsembleIntervals(i,2)),...
            yseq(IndxDixk),...
            'Color',ColorState(EnsembleInstances(i),:),'LineWidth',4);
        end
    end

    for i=1:Ntran
        % PLOT BALLS @ HEBBIAN SEQUENCE
        plot(EnsembleInstancesTimes(i)/fs/60,EnsembleInstances(i),'o',...
            'MarkerEdgeColor','k',...
            'MarkerFaceColor',ColorState(EnsembleInstances(i),:),...
            'MarkerSize',round(10*CellatEns(EnsembleInstances(i))+5));
        drawnow;
    end
    hold off;
    axis([Axis_details.XLim,0,max(EnsembleInstances)+1])
    EnsmbleAxis=gca;
    linkaxes([Axis_details,EnsmbleAxis],'x')
end