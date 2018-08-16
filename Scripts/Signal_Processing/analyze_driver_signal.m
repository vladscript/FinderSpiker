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
function [Dfix,XDfix,Xestfix,LambdasFix,IndexesFix,Features]=analyze_driver_signal(D,FR,XDupdate,Xest,varargin)
% check if Activate Driver Analysis
if isempty(varargin)
    check=1;
else
    check=0;
    disp('Only Checkin Driver Amplitudes')
end
% Make size of Cells x Frames:
D=makecellsxframes(D);
XDupdate=makecellsxframes(XDupdate);
Xest=makecellsxframes(Xest);
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
    disp(c);
    xd=XDupdate(c,:);
    xe=Xest(c,:);
    noisex=xd-xe;
    % Initial Maximum Driver->Distortion ##########################
    [~,framax]=max(d);
    if and(framax==1 && ~isempty(d(d~=0)) ,check)
        nextframe=framax+1;
        while and(nextframe<numel(xe),and(xe(nextframe)>0,x_sparse(nextframe)>0.5e-3))
            % x_sparse(nextframe)=0;
            % d(nextframe)=0;
            nextframe=nextframe+1;
        end
      
% Fix Initial Fast Bleaching
%         Apeaks=findpeaks(xe(1:nextframe));
%         if isempty(Apeaks)
%             XDfix(c,1:nextframe)=XDfix(c,1:nextframe)-x_sparsePRE(1:nextframe);
%         else
        if nextframe>3
            xexp=fit([1:nextframe]',xd(1:nextframe)','exp1');
            xdecay=xexp(1:nextframe);
            xd(1:nextframe)=xd(1:nextframe)-xdecay';
            XDfix(c,:)=xd;
            disp('\___ Initial decay ~~>')
        end
        
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
    % Thd=abs(min(d));
    SamplesDelete=1:numel(x_sparse);
    [Apeaks,Npeaks]=findpeaks(x_sparse);
    SaveSamples=[];
    % dx_sparse=diff(x_sparse);
    for n=1:length(Npeaks)
        if Apeaks(n)>=std(noisex)
            % Before the peak $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            auxN=0;
            isPeakStop=1;
            while and(and(x_sparse(Npeaks(n)-auxN)>0,Npeaks(n)-auxN>0),isPeakStop)
                if n>1
                    if and(ismember(Npeaks(n)-auxN,Npeaks(1:n-1)),x_sparse(Npeaks(n)-auxN)<std(noisex))
                        isPeakStop=0;
                    else
                        SaveSamples=[SaveSamples,(Npeaks(n)-auxN)];
                    end
                else
                    SaveSamples=[SaveSamples,(Npeaks(n)-auxN)];
                end
                auxN=auxN+1;
            end
            % After the peak $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
            auxN=1;
            isPeakStop=1;
            while and(and(x_sparse(Npeaks(n)+auxN)>0,Npeaks(n)+auxN<numel(x_sparse)),isPeakStop)
                if n<length(Npeaks)
                    if and(ismember(Npeaks(n)+auxN,Npeaks(n+1:end)),x_sparse(Npeaks(n)+auxN)<std(noisex))
                        isPeakStop=0;
                    else
                        SaveSamples=[SaveSamples,(Npeaks(n)+auxN)];
                    end
                else
                    SaveSamples=[SaveSamples,(Npeaks(n)+auxN)];
                end
                auxN=auxN+1;
                %lalalalala
                % SaveSamples=[SaveSamples,(Npeaks(n)+auxN)];
                % auxN=auxN+1;
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
    % Just (+) Drivers
    d(d<0)=0;
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
