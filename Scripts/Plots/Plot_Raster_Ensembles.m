%% Function to plot raster
% input
%   R:   Raster Matrix  Cells x Frames
%   varargin default:
%           fs=1            Sampling Frequency
%           stepy=2;        Step between rows Yticks
%           indexes=1:C;    Sort
%           ColocateIndx=[];Indexes of Colocalized Neurons
%           Neurons_Colocalized=[];
%   Neurons:Colocalized
% Output
% figure of the raster in MINUTES
% Always creates a new figure
function Plot_Raster_Ensembles(R,varargin)
%% Setup
TransientHight=0.75;
N=numel(varargin);
[C,~]=size(R);
switch N
    case 0
        fs=1; % Discrete Time
        stepy=2;
        indexes=1:C;
        ColocateIndx=[];% Indexes of Colocalized Neurons
        Neurons_Colocalized=[];
    case 1
        fs=varargin{1};
        stepy=2;
        indexes=1:C;
        ColocateIndx=[];% Indexes of Colocalized Neurons
        Neurons_Colocalized=[];
    case 2
        fs=varargin{1};
        stepy=varargin{2};
        indexes=1:C;
        ColocateIndx=[];% Indexes of Colocalized Neurons
        Neurons_Colocalized=[];
    case 3
        fs=varargin{1};
        stepy=varargin{2};
        indexes=varargin{3};
        ColocateIndx=[];% Indexes of Colocalized Neurons
        Neurons_Colocalized=[];
    case 4 
        fs=varargin{1};
        stepy=varargin{2};
        indexes=varargin{3};
        ColocateIndx=varargin{4}; % Indexes of Colocalized Neurons
        Neurons_Colocalized=find(ColocateIndx);
    otherwise
end

figure
ax1=subplot(3,1,[1,2]);
hold on;

if fs==1
    ts=1/fs*60;
else
    ts=1/fs;
end
% Re-Sort:
R=R(indexes,:);
for i=1:C
    j=indexes(i);
    
    ypositon=[j-TransientHight/2,TransientHight];
    activeframes=find(R(j,:));
    activeframes=[activeframes,0];
    % never and active frame is going to be Zero
    if ~isempty(activeframes)
        nf=1;
        % xposition(1)=activeframes(1)-0.5;
        xposition(1)=activeframes(1);
        xposition(2)=1;
        
        while nf<numel(activeframes)
            nx=nf;
            while activeframes(nx+1)==activeframes(nx)+1
                xposition(2)=xposition(2)+1;
                if nx+1==numel(activeframes)
                    % activeframes=[activeframes,0];
                    % never and active frame is going to be Zero
                    % Stop sloop;
                else
                    nx=nx+1;
                end
            end
            % Create Rectangle *************************
            xposs=ts*xposition/60;
            rectangle('Position',[xposs(1),ypositon(1),...
                    xposs(2),ypositon(2)],'Curvature',0.2,...
                    'EdgeColor','k',...
                    'FaceColor','k');
            fprintf('*')
            % Restart xposition values
            xposition(1)=activeframes(nx+1)-0.5;
            xposition(2)=1;
            nf=nx+1;
        end
        fprintf('\n')
    else
        fprintf('\n')
    end
    % If there is Colocalized Cells
    if ismember(j,Neurons_Colocalized)
            plot(ts*find(R(i,:))/60,i*R(i,R(i,:)>0)+0.6,'Marker','square',...
            'LineStyle','none',...
            'MarkerSize',2,...
            'MarkerEdgeColor',[1,0.6,0.78],...
            'MarkerFaceColor',[1,0.6,0.78]); hold on;
    end
end
axis([0,ts*(length(R))/60,1-TransientHight/2,C+TransientHight/2])
ylabel('Neural Activity')
set(gca,'XTick',[])
set(gca,'Box','off')
% Sorted Neurons according to Ensembles
yticks=indexes';
set(gca,'YTick',1:stepy:C)
set(gca,'YTickLabel',yticks(1:stepy:C))
set(gca,'TickLength',[0,0])

ax2=subplot(3,1,3);
plot(ts*(1:length(sum(R)))/60,sum(R),'k','LineWidth',1.1)
if ~isempty(ColocateIndx)
    hold on;
    plot(ts*(0:length(sum(R))-1)/60,sum(R(Neurons_Colocalized,:)),'Color',[1,0.6,0.78],'LineWidth',1.1)
end
hold off;

axis([0,ts*(length(R)-1)/60,0,max(sum(R))+1])
% grid on
ylabel('CAG')
if fs==1
    xlabel('Frames')
else
    xlabel('Minutes')
    set(gca,'XTick',0:1:ts*(length(R)/60))
end


set(gca,'Box','off')

lisp=linspace(round(max(sum(R))/2),max(sum(R)),2);
if lisp(1)~=lisp(2)
    set(gca,'YTick',floor(linspace(round(max(sum(R))/2),max(sum(R)),2)))
end
linkaxes([ax1,ax2],'x')
end