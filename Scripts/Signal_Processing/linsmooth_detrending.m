% Detrending by Segments Linear & Smoothing Methods
% Input
%   X: signals of Ns x Frames
% Output
%   XD: Detrended Signals
% function XD=linsmooth_detrending(X)
[Ns,Frames]=size(X); % Ns-Signals of Frames
% Initialize Outputs:
XD=zeros(Ns,Frames);
%% MAIN LOOP
for i=1:Ns 
    x=X(i,:); % Read Single Fluorescence Signal
    disp(['Detrending ... [Signal: ',num2str(i),'/',num2str(Ns),']']);
    y=smooth(x,Frames,'rloess');    % Trend Component (Bleaching)
    xd1=x-y';                       % 1st Detrending
    disp('------------------RLOESS smoothing OK')
    % PDF of Detrended Signal:
    [pxd1,binxd1]=ksdensity(xd1, linspace(min(xd1),max(xd1),100));
    % Modes of PDF of Detrended Signal:
    [~,Bin_Mode]=findpeaks(pxd1,binxd1,'SortStr','descend');
    % De-Offset:
    xd1=xd1-Bin_Mode(1); % Deoffset
    %% FIXING detrending **********************************************
     [cl,~,mucl] = polyfit(1:Frames,x,1);
     y_lin=polyval(cl,1:Frames,[],mucl); % Linear Trending
     xdlin=x-y_lin;         % 2nd DeTrending
     y_err=y_lin-y';        % Comparing Trendings
     [~,Bf]=findpeaks(-abs(y_err));
     Nintervals=numel(Bf)+1;
     a=1;       % Interval a:b
     ydet=[];  % Detrended Output Initalized
     for ni=1:Nintervals
             if ni<Nintervals
                 b=Bf(ni);
             else
                 b=Frames;
             end
             Frames2Test=a:b;
             a=b+1;
             % pdf of segment of the detrended signals
             [pxd1,binxd1]=ksdensity(xd1(Frames2Test));
             [pxdlin,binxdlin]=ksdensity(xdlin(Frames2Test));
             % Modes Positions of the pdf's
             [~,Bd]=findpeaks(pxd1,binxd1);
             [~,Bl]=findpeaks(pxdlin,binxdlin);
             if numel(Bd)==1 && numel(Bl)       % Case: One-Mode
                 disp('     Fixing detrending Linear-Smoothing:A')
                 [~,indxSEGMENT]=min(abs([Bd,Bl]));  % Get Closet to Zero
                 if indxSEGMENT==1                   % Chose Detrending Method
                        ydet=[ydet,y(Frames2Test)']; % Smoohting
                 else
                        ydet=[ydet,y_lin(Frames2Test)];    % Linear Detrendind
                 end
             else                               % Case: Multi Modes
                disp('     Fixing detrending Linear-Smoothing:B')
                % Compare Modes od Detrended Signal Closest to Zero
                if numel(Bd)>1
                    [~,indxD]=min(abs(Bd));
                    IndxModesD=setdiff(1:length(Bd),indxD);
                    ModesD=Bd(IndxModesD);      % Farthest Mode(s)
                else
                    ModesD=Bd;
                end
                if numel(Bl)>1
                    [~,indxL]=min(abs(Bl));
                    IndxModesL=setdiff(1:length(Bl),indxL);
                    ModesLin=Bl(IndxModesL);    % Farthest Mode(s)
                else
                    ModesLin=Bd;
                end 
                % If Negative Modes at Detrending Components:
                MD=0; ML=0;
                if isempty(ModesD(ModesD<0))
                    disp('No Negative Modes at smoothing Detrending')
                    MD=1;
                end
                if isempty(ModesLin(ModesLin<0))
                     disp('No Negative Modes at Linear Detrending')
                     ML=1;
                end
                 % Evaluate Different Cases 1 means acceptable
                 if MD==0 && ML==0 % Both Negative
                     disp('------------Warning A->')
                     DD=min(abs(Bd));
                     LL=min(abs(Bl));
                     [~,indxSEGMENT]=min([DD,LL]); % Closest Mode to Zero
                     if indxSEGMENT==1              % Chose Detrending Method
                            ydet=[ydet,y(Frames2Test)'];      % Smoohting
                     else
                            ydet=[ydet,y_lin(Frames2Test)];    % Linear Detrendind
                     end
                 elseif MD==0 && ML==1
                     disp('------------Warning B->')
                     ydet=[ydet,y_lin(Frames2Test)];    % Linear Detrendind
                 elseif MD==1 && ML==0
                     disp('------------Warning C->') 
                     ydet=[ydet,y(Frames2Test)'];      % Smoohting
                 elseif MD==1 && ML==1
                     disp('------------Warning D->')
                     DD=min(abs(Bd));
                     LL=min(abs(Bl));
                     [~,indxSEGMENT]=min([DD,LL]); % Closest Mode to Zero
                     if indxSEGMENT==1              % Chose Detrending Method
                            ydet=[ydet,y(Frames2Test)'];      % Smoohting
                     else
                            ydet=[ydet,y_lin(Frames2Test)];    % Linear Detrendind
                     end
                 end
             end
     end
     xd1=x-ydet;
     y=ydet;         
    % Preview Results
    subplot(211)
    plot(x); hold on;
    plot(y); hold off;
    axis tight; grid on;
    subplot(212)
    plot(xd1); hold on;
    plot([0,Frames],[0,0],'-.r');  hold off;
    axis tight; grid on;
    pause                 
    %% OUTPUT
    XD(i,:)=xd1;
end % function