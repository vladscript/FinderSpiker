% for i=1:Ns
%     x=X(i,:);
%     fc=0.015;
%     Norder=7;
%     [b,a]=butter(Norder,fc/(fs/2),'high');
%     xf=filtfilt(b,a,x); %detrended
%     y=smooth(xf,Frames,'rloess');
%     subplot(2,1,1)
%     plot(x);
%     subplot(2,1,2)
%     plot(xf,'k'); hold on;
%     plot(y,'k'); hold off;
%     pause;
% end
%% stuff
for i=1:Ns
    x=X(i,:);
    xd=XDupdate(i,:);
    xe=Xest(i,:);
    y=x-xd;
    noisex=xd-xe;
    subplot(2,1,1)
    plot(x); hold on;
    plot(y); hold off;
    axis tight; grid on;
    subplot(2,1,2)
    plot(xd,'k'); hold on;
    plot(xe,'m');
    plot([0,numel(x)],-[std(noisex),std(noisex)],'-.r');
    hold off;
    axis tight; grid on;
    pause;
end