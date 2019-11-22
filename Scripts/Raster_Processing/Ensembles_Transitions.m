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
if ~isempty(R)
    EnsemblesIndex=unique(labels_frames);
    for ne=1:numel(EnsemblesIndex)
        CellatEns(ne)= numel(find(sum(R(:,signif_frames( labels_frames==EnsemblesIndex(ne) )),2)))/size(R,1);
    end
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
    
    fig_Ensembles_Transitions=figure('Position', [415 172 560 122],...
        'Name','Hebbian Pathways');
    plot(EnsembleInstancesTimes/fs/60,EnsembleInstances,'k','LineWidth',0.5); hold on
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
        if ~isempty(EnsembleIntervals)
            % Get SLOPES
            tlin=[EnsembleIntervals(i,1):EnsembleIntervals(i,2)]/fs/60;
            if i==1
                % PRE SLOPE
                y2pre=EnsembleInstances(i+1);
                y1pre=EnsembleInstances(i);
                x2pre=EnsembleInstancesTimes(i+1)/fs/60;
                x1pre=EnsembleInstancesTimes(i)/fs/60;
                % POST SLOPE
                y2pos=y2pre;
                y1pos=y1pre;
                x2pos=x2pre;
                x1pos=x1pre;
            elseif i==Ntran
                % PRE SLOPE
                y2pre=EnsembleInstances(i);
                y1pre=EnsembleInstances(i-1);
                x2pre=EnsembleInstancesTimes(i)/fs/60;
                x1pre=EnsembleInstancesTimes(i-1)/fs/60;
                % POST SLOPE
                y2pos=y2pre;
                y1pos=y1pre;
                x2pos=x2pre;
                x1pos=x1pre;
            else
                % PRE-SLOPE
                y2pre=EnsembleInstances(i);
                y1pre=EnsembleInstances(i-1);
                x2pre=EnsembleInstancesTimes(i)/fs/60;
                x1pre=EnsembleInstancesTimes(i-1)/fs/60;
                % POST-SLOPE
                y2pos=EnsembleInstances(i+1);
                y1pos=EnsembleInstances(i);
                x2pos=EnsembleInstancesTimes(i+1)/fs/60;
                x1pos=EnsembleInstancesTimes(i)/fs/60;
            end
            % SLOPES
            pre_m=(y2pre-y1pre)/(x2pre-x1pre);
            pos_m=(y2pos-y1pos)/(x2pos-x1pos);
            % Y- INTERCEPT
            ytest=EnsembleInstances(i);
            xtest=EnsembleInstancesTimes(i)/fs/60;
            pre_b=ytest-pre_m*xtest;
            pos_b=ytest-pos_m*xtest;
            % LINE
            ypre=pre_m*tlin(tlin<=xtest)+pre_b;
            ypos=pos_m*tlin(tlin>xtest)+pos_b;
            % Plot Intervals
            plot(tlin(tlin<=xtest),ypre,'Color',ColorState(EnsembleInstances(i),:),'LineWidth',2)
            plot(tlin(tlin>xtest),ypos,'Color',ColorState(EnsembleInstances(i),:),'LineWidth',2)
            % PLOT @ HEBBIAN SEQUENCE
            plot(EnsembleInstancesTimes(i)/fs/60,EnsembleInstances(i),'o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',ColorState(EnsembleInstances(i),:),...
                'MarkerSize',round(10*CellatEns(EnsembleInstances(i))+5));
        else
            % PLOT @ HEBBIAN SEQUENCE
            plot(EnsembleInstancesTimes(i)/fs/60,EnsembleInstances(i),'o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor',ColorState(EnsembleInstances(i),:),...
                'MarkerSize',5);
        end
        %
    end
    hold off;
    axis([Axis_details.XLim,0.5,max(EnsembleInstances)+0.5])
    EnsmbleAxis=gca;
    linkaxes([Axis_details,EnsmbleAxis],'x')
end