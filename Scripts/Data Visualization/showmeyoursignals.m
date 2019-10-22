% Function to plot Calcium Transiens clean and raw
% in a sliding window for a simultaneous video (other script)
function showmeyoursignals(X_SPARSE,XD,Raster,SlideWin,NameOutSignals)

%%Setup
% if isempty(varargin{1})
%     SlideWin=500;
% else
%     SlideWin=varargin{1};
% end
disp('>>Started to Visualize:')
[Ns,Frames]=size(X_SPARSE);
% OUTPUT:
v = VideoWriter([NameOutSignals,'.avi']);
% v.FrameRate=fs; % Not really necessary
%     gray=0.8,0.8,0.8
% OFFSET SIGNALS:
offsetplot=0;
deltaoffset=0.2;
for i=1:Ns
    % Get Clean Signals & Normalize
    CellOffset(i)=offsetplot;
    CMax=max(X_SPARSE(i,:));
    xs=X_SPARSE(i,:)/CMax;
    xs(xs<0)=0;
    % Detrended Signal
    X_SPARSE(i,:)=xs;
    xd=XD(i,:)/CMax;
    % Update
    XD(i,:)=xd+offsetplot;
    X_SPARSE(i,:)=xs+offsetplot;
    offsetplot=offsetplot+max(xs)+deltaoffset;
end
signalswindow=figure; hold on;
LineTime=line([0,0],[0,offsetplot]);
LineTime.Color='k';

signalswindow.Children.YTick=[];
plot(XD','color',[0.45,0.45,0.45],'LineWidth',3); 
plot(X_SPARSE','color',[0,0.5,0],'LineWidth',1)
signalswindow.Children.XLim=[-SlideWin,0];
signalswindow.Children.YLim=[floor(min(XD(:))),ceil(max(XD(:)))];

open(v);
% offsetplot=0;
aux_count=1;
for f=1:Frames
    % Plot Time Line
    signalswindow.Children.XLim=[-SlideWin+f,0+f];
    LineTime.XData=[f,f];
    % Plot Spikes #################################
    % f=226
    for c=1:Ns
        if Raster(c,f)>0
            % Plot Line
            line([f,f],[CellOffset(c)-0.1,CellOffset(c)+0.1],'Color','k')
            % Plot Rectangle
        end
        % offsetplot=offsetplot+1;
        aux_count=aux_count+1;
    end
    
    Fmovie=getframe(signalswindow);
    writeVideo(v,Fmovie);
    fprintf('Recording %3.2f%% Video \n',100*f/Frames);
%     pause;
end
% hold off;  
close(v);
close(signalswindow);
disp('>>Done.')