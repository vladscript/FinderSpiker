% Funtion to Detrend Fluorescence Signals
% Input
%   X: signals of Ns x Frames
%   fs: sampling frequency
% Output
%   XD: Detrended Signals
function XD=only_detrending(X,fs)
[Ns,Frames]=size(X); % N-Signals of length Frames
ts=1/fs;                        % Sampling Period
t=linspace(0,Frames*ts,Frames); % Time Vector
% Initialize Outputs:
XD=zeros(Ns,Frames);
for i=1:Ns % Main LOOP
    x=X(i,:);
    % DETECTION OF FLUO DUE TO BASAL FIRING
    % pdf bimode detection 'fitgmdist'
    % if each mode belongs to different set of samples    
     % RLOESS   SMoothing Detrending
    disp(['Detrending ... [Signal: ',num2str(i),'/',num2str(Ns),']']);
    y=smooth(x,Frames,'rloess');                    % 1st Trend Component
    xd1=x-y';                                       % 1st Detrending
    
    % RLOWESS   SMoothing Detrending
    y=smooth(xd1,Frames,'rloess');                  % 2st Trend Component
    xd2= xd1-y';                                    % 2nd Detrending
    disp('Detrending ... {OK}')
    pxd1=ksdensity(xd1);
    Apeak=findpeaks(pxd1,'SortStr','descend','NPeaks',1);
    pxd2=ksdensity(xd2);
    Bpeak=findpeaks(pxd2,'SortStr','descend','NPeaks',1);
    if Apeak<Bpeak
        disp('>>>>>>>>>>>>>>>>>>>>>>>> Detrended OK .... ')
    else
        disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>Fix detrending .... ')
        % PDF of xd1 (approx) **************************************************
        xd1pts=linspace(min(xd1),max(xd1),100);
        [pxd1,binxd1]=ksdensity(xd1,xd1pts);
        [Pxd1_val,BinXd1,Wd1,~]=findpeaks(pxd1,binxd1);
        [~,NumMode]=max(Pxd1_val); % Biggest mode
        Th_Xd1=BinXd1(NumMode)+Wd1(NumMode)/2;
    
        % 2nd Order detrending of XD1 *********************************************
        samples_underTh=find(xd1<Th_Xd1);
        [cl,~,mucl] = polyfit(t(samples_underTh),xd1(samples_underTh),2);
        xd_poly=polyval(cl,t,[],mucl);
    %     xd_poly=smooth(xd1,Frames,'sgolay',2);    
        xd2=xd1-xd_poly;    % DETRENDED!  
    end
    XD(i,:)=xd2;         % SAVE DETRENDED * * * * ** $ $ $ $ OUTPUT
end