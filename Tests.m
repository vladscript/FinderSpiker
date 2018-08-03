% Derivative detrending
[C,F]=size(X);
for c=1:C
    x=X(c,:);
    
    plot(x); 
%     plot(y);
%     hold off;
    axis tight; grid on;
    
    pause
end

% Filtering

% AL New_Index are reffered to Original Data: SIGNALS, DETSIGNALS, RASTEROK XY
% New_Index_Active=New_Index(TotalActiveNeurons); % Sorted Indexes of Clean Rasters
% [RASTER_Selected_Clean,XY_selected,R_Condition,Onsets,New_Index_Sel]= Select_Raster_for_NN(fs,Raster_Condition,XY,Names_Conditions,Experiment);
% 
% [XY_merged,ColocateIndx]=get_merged_coordinates(Experiment,XY_clean,r);
% A2A=find(ColocateIndx); % 
% Indexes=New_Index(A2A);
% 
% %%CORRELATION***********************************************
% simindex=corrcoef(RASTER_Selected_Clean');
% cluster_index=clusterModularity(simindex,1000);
% plotClusterRaster(RASTER_Selected_Clean,cluster_index,1);

%%% Plot Some PDFs to monitor errors (!)
%             figure; 
%             subplot(211)
%             if and(~isempty(AcceptedINDX),~isempty(RejectedINDX))
%                 [pA,binA]=ksdensity(SNRbyWT(AcceptedINDX),linspace(min(SNRbyWT(AcceptedINDX)),max(SNRbyWT(AcceptedINDX)),100));
%                 [pR,binR]=ksdensity(SNRbyWT(RejectedINDX),linspace(min(SNRbyWT(RejectedINDX)),max(SNRbyWT(RejectedINDX)),100));
%                 plot(binR,pR,'.r','LineWidth',2)
%                 hold on;
%                 plot(binA,pA,'*b','LineWidth',2)
%                 hold off;
%                 legend('Rejected','Accepted','Location','northwest')
%             else
%                 [pO,binO]=ksdensity(SNRbyWT,linspace(min(SNRbyWT),max(SNRbyWT),100));
%                 plot(binO,pO,'.m','LineWidth',1)
%             end
%             axis tight; grid on;
%             title('SNR pdf')
%             subplot(212)
%             if and(~isempty(AcceptedINDX),~isempty(RejectedINDX))
%                 [pA,binA]=ksdensity(SkewSignal(AcceptedINDX),linspace(min(SkewSignal(AcceptedINDX)),max(SkewSignal(AcceptedINDX)),100));
%                 [pR,binR]=ksdensity(SkewSignal(RejectedINDX),linspace(min(SkewSignal(RejectedINDX)),max(SkewSignal(RejectedINDX)),100));
%                 plot(binR,pR,'.r','LineWidth',2)
%                 hold on;
%                 plot(binA,pA,'*b','LineWidth',2)
%                 hold off;
%                 legend('Rejected','Accepted','Location','northwest')
%             else
%                 [pO,binO]=ksdensity(SkewSignal,linspace(min(SkewSignal),max(SkewSignal),100));
%                 plot(binO,pO,'.m','LineWidth',1)
%             end
%             axis tight; grid on;
%             title('Skewness pdf')
% figure
% % for i=1:Ns 
%     x=X(i,:); % Read Single Fluorescence Signal
%     % if each mode belongs to different set of samples    
%      % RLOESS   SMoothing Detrending
%     disp(['Detrending ... [Signal: ',num2str(i),'/',num2str(Ns),']']);
%     y=smooth(x,Frames,'rloess');                    % Trend Component (Bleaching)
%     xd1=x-y';
%     %% Test Detrending before Distortion
%     [~,PeaksY]=findpeaks(y);
%     [~,ValleY]=findpeaks(-y);
%     [PeaksY,ValleY]
%     [px,binx]=ksdensity(x, linspace(min(x),max(x),100));
%     [pxd1,binxd1]=ksdensity(xd1, linspace(min(xd1),max(xd1),100));
%     subplot(2,3,[1,2])
%     plot(x); hold on;
%     plot(y); hold off;
%     subplot(2,3,[4,5])
%     plot(xd1)
%     axis tight; grid on;
%     subplot(2,3,3)
%     plot(px,binx);
%     grid on;
%     subplot(2,3,6)
%     plot(pxd1,binxd1);
%     % axis tight;
%     grid on;
%     pause;
% end
% %% PLOT SIGNALS
% 
% for i=1:188
%     disp(i)
%     subplot(3,1,1)
%     plot(X(i,:));
%     axis tight; grid on;
%     subplot(3,1,2)
%     plot(XD(i,:),'b'); hold on;
%     plot(XDupdate(i,:),'r'); hold off;
%     axis tight; grid on;
%     subplot(3,1,3)
%     if i==okINDX(i)
%         plot(XDfix(i,:),'g');
%         axis tight; grid on;
%     end
%     pause
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
