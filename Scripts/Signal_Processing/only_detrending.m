% Funtion to Detrend Fluorescence Signals
% Input
%   X: signals of Ns x Frames
% Output
%   XD: Detrended Signals
function XD=only_detrending(X)
[Ns,Frames]=size(X); % Ns-Signals of Frames
% Initialize Outputs:
XD=zeros(Ns,Frames);
%% MAIN LOOP
for i=1:Ns 
    x=X(i,:); % Read Single Fluorescence Signal
    % if each mode belongs to different set of samples    
     % RLOESS   SMoothing Detrending
    disp(['Detrending ... [Signal: ',num2str(i),'/',num2str(Ns),']']);
    y=smooth(x,Frames,'rloess');    % Trend Component (Bleaching*) *BEST*
    % y=smooth(x,Frames,'loess');    % Trend Component (Bleaching*)
    %% Test Detrending before Distortion
    % PDF of Raw Signal:
    [px,binx]=ksdensity(x, linspace(min(x),max(x),100));
    [PModesX,ModesX]=findpeaks(px,binx,'SortStr','descend');
    ModesXaux=sort(ModesX);
    if ModesXaux(1)==ModesX(1) % Check if biggest is least gray
        yaux=smoothbysegments(x,y,1); % More than 1 Curve
        xd1=x-yaux;   
        disp('------------------RLOESS smoothing OK [1]')
    else % If first mode ain't the biggest one
        yaux=smoothbysegments(x,y,0); % Any Curve
        xd1=x-yaux;   
        disp('------------------RLOESS smoothing OK [0]')
    end
    % Second Detrending Test
    [~,ylinv1]=lineardetbigmode(binx,ModesX,x);
    xlin=x-ylinv1;

    % PDF of Detrended Signal:
    [pxd1,binxd1]=ksdensity(xd1, linspace(min(xd1),max(xd1),100));
    % Modes of PDF of Detrended Signal:
    [~,Bin_Mode]=findpeaks(pxd1,binxd1,'SortStr','descend');
    % De-Offset:
    xd1=xd1-Bin_Mode(1); % Deoffset Biggest Mode
    % PDF of Detrended Signal RECALCUL
    [pxd1,binxd1]=ksdensity(xd1, linspace(min(xd1),max(xd1),100));
    % Modes of PDF of Detrended Signal:
    [Mode_Peak,Bin_Mode]=findpeaks(pxd1,binxd1,'SortStr','descend');
    %% CHECK FOR BASAL ACTIVITY
    if numel(ModesX)>1 && ModesX(1)==ModesXaux(1)
        NegatiVeDetrended=find(xlin<0);
        NegatiVeDetrendedSmooth=find(xd1<0);
        HiguerSamples=find(x>max(ModesX(1:2)));
    else
        disp('>> Ok Detrended by RLOESS smoothing');
    end
    %% MissDetrending Checking 2nd Part *
    if numel(Mode_Peak)>1
        % Considerations:
        % If there are Negative Modes-> Missdetrending Issues
        % Biggest Mode may be Negative but close to Zero, so ingore it:
        if sum(Bin_Mode(2:end)<0)==0
            disp('>> Possible Calcium Transients Detected')
        else
            disp('>> MissDetrending Flag [ ]')
            %% FIXING detrending 
             [cl,~,mucl] = polyfit(1:Frames,x,1);
             y_lin=polyval(cl,1:Frames,[],mucl); % Linear Trending
             xdlin=x-y_lin;         % Linear DeTrending
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
                    % Plot preview
                 
                 end
                 
             end
             xd1=x-ydet;
             y=ydet;
        end
    else 
        disp('ok...')

    end
    % Preview Results
    % CHECK FIGURE
    [pxd1,binxd1]=ksdensity(xd1, linspace(min(xd1),max(xd1),100));
    y=x-xd1;
    Sup=subplot(2,3,[1,2]);
    plot(x); hold on;
    plot(y); 
    axis tight; grid on;
    hold off;
    Sdown=subplot(2,3,[4,5]);
    plot(xd1)
    axis tight; grid on;
    Hup=subplot(2,3,3);
    plot(px,binx);
    grid on;
    Hdown=subplot(2,3,6);
    plot(pxd1,binxd1);
    axis tight;
    grid on;
    linkaxes([Sup,Sdown],'x');
    linkaxes([Sup,Hup],'y');
    linkaxes([Sdown,Hdown],'y');
    disp('anfubvtqpsm')
    
    %% OUTPUT
    XD(i,:)=xd1;
end % loop
end % Main Function
%% Nested Function
function [samplestodetrend,ylinv1]=lineardetbigmode(binx,ModesX,x)
    % Get Samples belonging from least sample to the biggest mode
    Frames=length(x);
    ModePos=find(binx==ModesX(1));
    % dpx=diff(px);
    samplestodetrend= intersect( find(x>min(x)),find(x<binx(ModePos)) );
    [cl,~,mucl] = polyfit(samplestodetrend,x(samplestodetrend),1);
    ylinv1=polyval(cl,1:Frames,[],mucl); % Linear Trending
end

function yss = smoothbysegments(x,y,NC)
    [~,PeaksY]=findpeaks(y);
    [~,ValleY]=findpeaks(-y);
    Curves=[makerowvector(PeaksY),makerowvector(ValleY)];
    if numel(Curves)>NC
        disp('   -^-v- Curvy Trending )))')
        Curves=sort(Curves);
        frameA=1;
        yss=zeros(size(x));
        for n=1:numel(Curves)+1
            if n<=numel(Curves)
                frameB=Curves(n);
            else
                frameB=numel(x);
            end
            xseg=x(frameA:frameB);
            yseg=smooth(xseg,numel(xseg),'rloess'); 
            % (RE)-CURSIVENESS:############ NOT YET
            % yseg=smoothbysegments(xseg,yseg);
            yss(frameA:frameB)=yseg;
            frameA=frameB;
        end
        
    else
        yss=y';
    end
end