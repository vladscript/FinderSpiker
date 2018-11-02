%% Function to Plot Ensembles **************************************
% Input
%   EnsembleName        Condition Names
%   Ensembled_Raster    Raster
%   Ensembled_Labels    Labels in time
%   Ensemble_Threshold  Significane Threshold
%   UniRMutiE           Unique Raster Multi Conditions Indicator
%   ColorState          Color Map 
%   fs                  Sampling Frequency
% Output
%   Figure of Ensembles by Condition
function Plot_Ensembles(EnsembleName,Ensembled_Raster,Ensembled_Labels,...
        Ensemble_Threshold,UniRMutiE,ColorState,fs)
%% Setup
NCplot=numel(EnsembleName);
AuxC=0;
%% Plot Whole Experiment
%% Conditions Loop
for c=1:NCplot
    R=Ensembled_Raster{c};      % Frames x Cels
    labels=Ensembled_Labels{c}; % Clustering Labels
    CoAc=sum(R,2);             % Coactivity Signal
    THR{c}=Ensemble_Threshold{c};
    sigframes=find(CoAc>=THR{c}); % Significative Frames
    NG=numel(unique(labels));           % N ensembles @ condition 
    LENGHTRASTER{c}=length(Ensembled_Raster{c});          % Raster Length
%     ActNeu=find(sum(R,1)>0);   % Active Neurons in each Raster: ENS & NON ENS
    if NG>0     % IF THERE'RE ENSEMBLES
        % Re-Sorting in each Condition:
        [OrderOneCondition,~]=OrderClusters(labels,sigframes,R,NG);
        % SORTING CONDITION
        Index_Ensemble=OrderOneCondition;          % Neurons Label Raster
        % Plotting |
        Plot_Raster_Ensembles(R,Index_Ensemble,5,fs);       % Sorted Raster
        EnsembleFig=gcf; EnsembleFig.Name=EnsembleName{c};
        if UniRMutiE
            % Included (deep) Purple:
            CS=ColorState;
        else
            CS = ColorState(AuxC+1:AuxC+NG,:);          
            CS=[CS;ColorState(end,:)];  % Plus (deep) Purple
        end
        Plot_State_Colors(labels,sigframes,CS,R,fs,CoAc,Index_Ensemble);        
        plot_CAG_threshold(THR,LENGHTRASTER,fs);
        % Cycles Reverberation Analysis: in waiting
        drawnow;
        
    else            % Non-Ensembles
        
        Plot_Raster_V(R,fs);
        EnsembleFig=gcf; EnsembleFig.Name=EnsembleName{c};
        drawnow;
    end
    % Update auxiliar variables
    AuxC=AuxC+NG;
end
