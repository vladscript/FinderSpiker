% Plot detrended signals
function Plot_detrended(XD)
[Ns,Frames]=size(XD);
figure; hold on;
offsetplot=0;
for i=1:Ns
    xd=XD(i,:);
    if sum(xd)==0
        xd(:)=0.5;
        disp('Empty Signal ... ')
        offsetplot=offsetplot+1;
    else
        C=max(xd)-min(xd);
        xd=(xd-min(xd))/C;
        plot(xd+offsetplot,'k')
        offsetplot=offsetplot+max(xd);
    end
end
axis([0,Frames,0,Ns])
grid on
hold off;