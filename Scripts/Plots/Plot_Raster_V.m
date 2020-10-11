%% Function to plot raster
% Input
%   R: Raster Matrix Cells x Frames
%   fs= sampling frequency
% Output
%   Figure of the raster in MINUTES
%   Always creates a new figure
function Plot_Raster_V(R,varargin)
if isempty(varargin)
    fs=1/60; % frames
else
    fs=varargin{1};
end
figure
ax1=subplot(3,1,[1,2]);
hold on;
% % Identify Cell and Frame Matrix Dimension
[C,F]=size(R);
% if C>F
%     R=R';
%     [C,~]=size(R);
% end

ts=1/fs;
for i=1:C
    plot(ts*find(R(i,:))/60,i*R(i,R(i,:)>0),'o',...
                'MarkerEdgeColor','k',...
                'MarkerFaceColor','k',...
                'MarkerSize',3); hold on;
end
axis([0,ts*(length(R)-1)/60,0,C+1])
ylabel('Neural Activity')
set(gca,'XTick',[])
set(gca,'Box','off')
if C>9
    yticks=unique(  round(linspace(1,C,5)) );
    set(gca,'YTick',yticks(2:end))
else
    yticks=1:C;
    set(gca,'YTick',yticks)
end
set(gca,'TickLength',[0,0])

ax2=subplot(3,1,3);
plot(ts*(0:length(sum(R))-1)/60,sum(R),'k','LineWidth',1.1)
axis([0,ts*(length(R)-1)/60,0,max(sum(R))+1])
% grid on
ylabel('CAG')
if isempty(varargin)
    ax2.XLabel.String='Frames';
else
    ax2.XLabel.String='Minutes';
end

set(gca,'Box','off')
% set(gca,'XTick',0:1:ts*(length(R)/60))
lisp=linspace(round(max(sum(R))/2),max(sum(R)),2);
if lisp(1)~=lisp(2)
    set(gca,'YTick',floor(linspace(round(max(sum(R))/2),max(sum(R)),2)))
end
linkaxes([ax1,ax2],'x')
end