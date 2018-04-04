% Function to plot Calcium Transiens clean and raw
% in a sliding window for a simultaneous video (other script)
function showmeyoursiglas(X_SPARSE,Raster)
[Ns,Frames]=size(X_SPARSE);
v = VideoWriter([pwd,'CaTransients.avi']);
    % v.FrameRate=fs; % Not really necessary
    
% OFFSET SIGNALS:
offsetplot=0.2;
for i=1:Ns
    % Get Clean Signals & Normalize
    xs=X_SPARSE(i,:);
    xs(xs<0)=0;
    xs=0.8*xs/max(xs)+offsetplot;
    X_SPARSE(i,:)=xs;
    offsetplot=offsetplot+1;
end
signalswindow=figure; hold on;
LineTime=line([1,1],[0,Ns+1]);
LineTime.Color='k';
axis([0,Frames,0,Ns])
signalswindow.Children.YTick=1:Ns;
plot(X_SPARSE','b')
open(v);
% offsetplot=0;
aux_count=1;
for f=1:Frames
    % Plot Time Line
    LineTime.XData=[f,f];
    % Plot Spikes
    for c=1:Ns
        if Raster(c,f)>0
            % Plot Line
            line([f,f],[c-1,c-1+0.2],'Color','k')
        end
        % offsetplot=offsetplot+1;
        aux_count=aux_count+1;
    end
    
    Fmovie=getframe(signalswindow);
    writeVideo(v,Fmovie);
    disp(['Recording Video ... ',num2str(100*aux_count/(Frames*Ns),5),'%'])
%     pause;
end
% hold off;    