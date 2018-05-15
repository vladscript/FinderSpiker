% Function to anlize Driver Function:
% Check if first samples are the biggest ones (mainly due to miss
% detrending, fast bleaching or decay of fluorescence
% Input
%   D:              Original Set of Row Vector of Driver Signals
%   FR:             Response Functions
%   XDupdate:       Detrendend Signals
%   SigmaNoise:     Standard Deviation of Noise
% Output
%   Dfix:           Fixed Set of Row Vector of Driver Signals
%   XDfix:          Detrended Signals Fixed
%   Xestfix:        Denoised Signals Fixed
%   LambdasFix:     Lambdas Recalculated
%   IndexesFix:     Indexes of Fixed Cells
%   Features:       [SkewSignal, SkewNoise, SRbyWT] Matrix of Fetures]
function [Dfix,XDfix,Xestfix,LambdasFix,IndexesFix,Features]=analyze_driver_signal(D,FR,XDupdate,Xest)
% Initial Decay due to firing or "fast bleaching"
[C,~]=size(D); % Length of okINDX
% Initialize Output
Xestfix=Xest;
XDfix=XDupdate;
Dfix=D;
LambdasFix=[];
IndexesFix=[];
fixindx=1;
Features=[];
for c=1:C
    % Get Single Signals ##########################################
    r=FR(c,:);  % response function
    d=D(c,:);   % driver function
    x_sparse=sparse_convolution(d,r); % Clean Signal
    % x_sparsePRE=x_sparse;
    xd=XDupdate(c,:);
    xe=Xest(c,:);
    noisex=xd-xe;
    % Initial Maximum Driver->Distortion ##########################
    [~,framax]=max(d);
    if framax==1 && sum(d)~=0
        nextframe=framax+1;
        while and(xe(nextframe)>0,x_sparse(nextframe)>0.5e-3)
            % x_sparse(nextframe)=0;
            % d(nextframe)=0;
            nextframe=nextframe+1;
        end
      
% Fix Initial Fast Bleaching
%         Apeaks=findpeaks(xe(1:nextframe));
%         if isempty(Apeaks)
%             XDfix(c,1:nextframe)=XDfix(c,1:nextframe)-x_sparsePRE(1:nextframe);
%         else
        xexp=fit([1:nextframe]',xd(1:nextframe)','exp1');
        xdecay=xexp(1:nextframe);
        xd(1:nextframe)=xd(1:nextframe)-xdecay';
        XDfix(c,:)=xd;
        disp('~~>')
%         end
        [xe,noisex]=mini_denoise(xd); % update and fix
        [d,lambdaD]=maxlambda_finder(xd,r);
        x_sparse=sparse_convolution(d,r);
        Xestfix(c,:)=xe;
        % Dfix(c,:)=d; % update and fix
        LambdasFix=[LambdasFix,lambdaD];
        IndexesFix=[IndexesFix,c];
        % Recalculate Signal Features:
        [Features(fixindx,1),Features(fixindx,2),Features(fixindx,3),~]=feature_extraction(xe,noisex);
        fixindx=fixindx+1;
        disp('Initial fast bleaching or decaying')
    else
        disp('Fluoresence Trace OK')
    end
    % Delete Small Changes below Noise #############################
    % Check Responses below Noise:
    Thd=abs(min(d));
    d(d<=Thd)=0;
    SamplesDelete=1:numel(x_sparse);
    [Apeaks,Npeaks]=findpeaks(x_sparse);
    SaveSamples=[];
    for n=1:length(Npeaks)
        if Apeaks(n)>=std(noisex)
            % Before from the peak
            auxN=0;
            while and(x_sparse(Npeaks(n)-auxN)>0,Npeaks(n)-auxN>0)
                SaveSamples=[SaveSamples,(Npeaks(n)-auxN)];
                auxN=auxN+1;
            end
            % After from the peak
            auxN=1;
            while and(x_sparse(Npeaks(n)+auxN)>0,Npeaks(n)+auxN<numel(x_sparse))
                SaveSamples=[SaveSamples,(Npeaks(n)+auxN)];
                auxN=auxN+1;
            end
        end
    end
    SaveSamples=unique(SaveSamples);
    if ~isempty(SaveSamples)
        SamplesDelete=setdiff([1:numel(x_sparse)],SaveSamples);
        d(SamplesDelete)=0;
    else
        d(SamplesDelete)=0;
    end
    
    
    Dfix(c,:)=d; % update and fix
%     %% CHECK STUFF
%     plot(xd,'b'); hold on;
%     plot(xe,'m'); hold on;
%     plot(x_sparse,'g','LineWidth',2);
%     plot([0,numel(x_sparse)],[std(noisex),std(noisex)],'-.r');
%     bar(d); hold off;
%     %pause
%     disp(c) 
end
