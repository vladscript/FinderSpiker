% Function to make funstional Network from Raster Matrix
% Input
%   R:                  Condition Raster: Cells x Frames
%   XY:                 Original Coordinates: [x,y]
%   r:                  radius: vector
%   NameCondition:      Name of the Condition: String array
% Output
%   A: Adjacency Matrix
%   IndxNeuron: Index of Original Coordinates
%   NetFeatures: Network features:
%       gamma:  p(k) links distribution
%       alpha:  c(k) clustering index index
%       r^2:    fit
%   BannedFrames:   Relative Deleted Frames
% MANUAL HEYS COMMANDS:
%   r:  Reset Original Data
%   l:  Low Threshold Link Frames
%   u:  Up Threshold Link Frames
function [A,IndxNeuron,NetFeatures,BannedFrames]=Raster2fNet(R,XY,r,NameCondition)
%% SETUP*************************
[TotalCells,~]=size(R);
IndxNeuron=1:TotalCells;
Active_Neurons=find(sum(R,2)>0);        % Get Only Active Ones
IndxNeuron=IndxNeuron(Active_Neurons);
XY_active=XY(Active_Neurons,:);
radius_activus=r(Active_Neurons);
circle_areas=radius_activus.*radius_activus*pi;

% Initialize Data:
frames2del=[];              % Banned Frames
BannedFrames=[];            % ALL Banned Frames
th_frames=1;                % Threshold of Links - Frames
Rx=R(Active_Neurons,:);     % Raster from Only Active Cells
CoAc=sum(Rx);               % CoActivity Signal
[AN,Frames]=size(Rx);        % N frames
% [A,Pk,Ck]=Get_Network(Rx);
%% Initialize Output from Complete Raster    ******************
[A,kdegree,pk,binK,ck,NetFeatures]=NetParameters(Rx,th_frames);
klinks=unique(kdegree);
pk_fit=NetFeatures(1,1).*klinks.^-(NetFeatures(1,2));
ck_fit=NetFeatures(2,1).*klinks.^-(NetFeatures(2,2));
%% Main Figure: Initialize ++++++++++++++++++++++++++++++++++++++
checknetwork=figure('numbertitle','off',...
            'name',NameCondition,'position',[46 42 800 400],...
            'keypressfcn',@manual_cmnd);

% Define colormap
CM=hot(20);         % Choose Colors
CM=CM(end:-1:1,:);  % Turn it up-down
colormap(CM);       % Make Color map
%% Initialize AXIS and PLOTS (with randomw stuff) #########################

ax1=subplot(3,4,[1,2,5,6]);             % RASTER---------------------------
plotraster=imagesc(Rx,'Parent',ax1);
plotraster.Parent.Title.String='RASTER';
% plotraster.Parent.XTick=[];
plotraster.Parent.Box='off';
set(plotraster,'ButtonDownFcn',@bann_columns)

ax2=subplot(3,4,[9,10]); % Coactivity Signal ------------------------------
plotcoac=plot(CoAc,'Parent',ax2);
plotcoac.Parent.YLabel.String='Co-Activity';
plotcoac.Parent.XLim=[1,Frames];
plotcoac.Parent.YLim=[min(CoAc),max(CoAc)];
plotcoac.Parent.XGrid='on';
plotcoac.Parent.YGrid='on';
set(plotcoac,'ButtonDownFcn',@bann_columns)

ax3=subplot(3,4,[3,4,7,8]); % fNetwork @ Original Tissue Coordinates-------
set(ax3,'NextPlot','replacechildren');
hold(ax3,'on')
% Locations
scatter(ax3, XY_active(:,1),XY_active(:,2),circle_areas);
% Links
for n=1:AN
    LinkN=find(A(n,:));
    for l=1:length(LinkN)
        plot(ax3,[XY_active(n,1),XY_active(LinkN(l),1)],[XY_active(n,2),XY_active(LinkN(l),2)],...
            'k','LineWidth',3*A(n,LinkN(l))/max(max(A)));
    end
end
hold(ax3,'off')
title(ax3,'fNETWORK');
% set(plotnetwork,'ButtonDownFcn',@remove_driver)

ax4=subplot(3,4,11); % p(k)------------------------------------------------
set(ax4,'NextPlot','replacechildren');
hold(ax4,'on')
plot(ax4,binK(pk>0),pk(pk>0),'ok','MarkerSize',5)
plot(ax4,klinks,pk_fit,'b','LineWidth',2)
grid(ax4,'on');
hold(ax4,'off')
axis(ax4,[min(klinks),max(klinks),0,1]);
title(ax4,['\gamma=',num2str(NetFeatures(1,2),2),'  R^2=',...
    num2str(NetFeatures(1,3),2)],'FontSize',8,'FontWeight','normal');
% plotpk=plot(kdegree,pk,'Parent',ax4);     % Fit
% plotpk.Parent.YLabel.String='p(k)';
% plotraster.Parent.Box='off';

ax5=subplot(3,4,12); % c(k)------------------------------------------------
set(ax5,'NextPlot','replacechildren');
hold(ax5,'on')
plot(ax5,klinks(ck>0),ck(ck>0),'ok','MarkerSize',5)
plot(ax5,klinks,ck_fit,'b','LineWidth',2)
grid(ax5,'on');
hold(ax5,'off')
axis(ax5,[min(klinks),max(klinks),0,1]);
title(ax5,['\alpha=',num2str(NetFeatures(2,2),2),'  R^2=',...
    num2str(NetFeatures(2,3),2)],'FontSize',8,'FontWeight','normal');
% plotck=plot(kdegree,ck,'Parent',ax5);     % Fit
% plotck.Parent.YLabel.String='c(k)';
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% WAIT FOR OUTPUT
pause
%% NESTED FUNTIONS

%% MANUAL CTRL
    function manual_cmnd(checknetwork,~,~)
        key=get(checknetwork,'CurrentKey');
        disp(key)
        if strcmp(key,'r')
            [TotalCells,~]=size(R);
            IndxNeuron=1:TotalCells;
            Active_Neurons=find(sum(R,2)>0);        % Get Only Active Ones
            IndxNeuron=IndxNeuron(Active_Neurons);
            XY_active=XY(Active_Neurons,:);
            radius_activus=r(Active_Neurons);
            circle_areas=radius_activus.*radius_activus*pi;
            % Initialize Data:
            frames2del=[];              % Banned Frames
            BannedFrames=[];            % ALL Banned Frames
            th_frames=1;                % Threshold of Links - Frames
            Rx=R(Active_Neurons,:);     % Raster from Only Active Cells
            CoAc=sum(Rx);               % CoActivity Signal
            [AN,Frames]=size(Rx);       % N frames
            % [A,Pk,Ck]=Get_Network(Rx);
            %% Initialize Output from Complete Raster    ******************
            [A,kdegree,pk,binK,ck,NetFeatures]=NetParameters(Rx,th_frames);
            klinks=unique(kdegree);
            pk_fit=NetFeatures(1,1).*klinks.^-(NetFeatures(1,2));
            ck_fit=NetFeatures(2,1).*klinks.^-(NetFeatures(2,2));
            Re_Plot_Data_Now;
            disp('Reseted Raster')
        elseif strcmp(key,'u')
            th_frames=th_frames+1;                % Threshold of Links - Frames
            [A,kdegree,pk,binK,ck,NetFeatures]=NetParameters(Rx,th_frames);
            klinks=unique(kdegree);
            pk_fit=NetFeatures(1,1).*klinks.^-(NetFeatures(1,2));
            ck_fit=NetFeatures(2,1).*klinks.^-(NetFeatures(2,2));
            Re_Plot_Data_Now;
            disp('Rised Frames Threshold ')
        elseif strcmp(key,'l')
            if th_frames>2
                th_frames=th_frames-1;                % Threshold of Links - Frames
            else
                th_frames=1;
            end
            [A,kdegree,pk,binK,ck,NetFeatures]=NetParameters(Rx,th_frames);
            klinks=unique(kdegree);
            pk_fit=NetFeatures(1,1).*klinks.^-(NetFeatures(1,2));
            ck_fit=NetFeatures(2,1).*klinks.^-(NetFeatures(2,2));
            Re_Plot_Data_Now;
            disp('Lowed Frames Threshold ')
        end
    end

%% Replot Data
    function Re_Plot_Data_Now()
        % RASTER---------------------------
        plotraster.CData=Rx;
        plotraster.Parent.YLim=[0,AN];
        plotraster.Parent.XLim=[1,Frames];
        % plotraster.Parent.XTick=[];
        plotraster.Parent.Box='off';
        
        % Coactivity Signal ------------------------------
        plotcoac.YData=CoAc;
        axis(plotcoac.Parent,'tight');
        grid(plotcoac.Parent,'on');

        % fNetwork @ Original Tissue Coordinates-------
        cla(ax3);
        hold(ax3,'on')
        % Locations
        scatter(ax3, XY_active(:,1),XY_active(:,2),circle_areas);
        % Links
        for n=1:AN
            LinkN=find(A(n,:));
            for l=1:length(LinkN)
                plot(ax3,[XY_active(n,1),XY_active(LinkN(l),1)],[XY_active(n,2),XY_active(LinkN(l),2)],...
                    'k','LineWidth',3*A(n,LinkN(l))/max(max(A)));
            end
        end
        hold(ax3,'off')
        title(ax3,['fNETWORK Active Cells: ',num2str(AN)]);

        % p(k)------------------------------------------------
        cla(ax4);
        hold(ax4,'on')
        plot(ax4,binK(pk>0),pk(pk>0),'ok','MarkerSize',5)
        plot(ax4,klinks,pk_fit,'b','LineWidth',2)
        grid(ax4,'on');
        hold(ax4,'off')
        axis(ax4,[min(klinks),max(klinks),0,1]);
        title(ax4,['\gamma=',num2str(NetFeatures(1,2),2),'  R^2=',...
            num2str(NetFeatures(1,3),2)],'FontSize',8,'FontWeight','normal');
        

        % c(k)------------------------------------------------
        cla(ax5);
        hold(ax5,'on')
        plot(ax5,klinks(ck>0),ck(ck>0),'ok','MarkerSize',5)
        plot(ax5,klinks,ck_fit,'b','LineWidth',2)
        grid(ax5,'on');
        hold(ax5,'off')
        axis(ax5,[min(klinks),max(klinks),0,1]);
        title(ax5,['\alpha=',num2str(NetFeatures(2,2),2),'  R^2=',...
            num2str(NetFeatures(2,3),2)],'FontSize',8,'FontWeight','normal');
        drawnow;
        % +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end

%% Update Data: Run After Banning Frames
    function UpdateData()
        Rx=Rx(:,setdiff(1:Frames,[frames2del(1):frames2del(2)]));      % Bann Frames
        Active_Neurons=find(sum(Rx,2)>0);           % Get Only Active Cells:
        IndxNeuron=IndxNeuron(Active_Neurons);      % Original Neuron ID
        XY_active=XY(Active_Neurons,:);             % Active Coordinates
        radius_activus=r(Active_Neurons);           % Active Radius
        circle_areas=radius_activus.*radius_activus*pi;
        % New Raster
        Rx=Rx(Active_Neurons,:);     % Raster from Only Active Cells
        CoAc=sum(Rx);               % CoActivity Signal
        [AN,Frames]=size(Rx);        % N frames
        % [A,Pk,Ck]=Get_Network(Rx);
        %% Initialize Output from Complete Raster    ******************
        [A,kdegree,pk,binK,ck,NetFeatures]=NetParameters(Rx,th_frames);
        klinks=unique(kdegree);
        pk_fit=NetFeatures(1,1).*klinks.^-(NetFeatures(1,2));
        ck_fit=NetFeatures(2,1).*klinks.^-(NetFeatures(2,2));
    end
%% Bann Column Vectors
    function bann_columns(~,~)
         % rect: [xmin ymin width height]
        rect=getrect; rect=round(rect);
        frames2del=[rect(1),rect(1)+rect(3)];
        % Limit Window to Delete
        if frames2del(1)<1
            frames2del(1)=1;
        end
        if frames2del(2)>Frames
            frames2del(2)=Frames;
        end
        UpdateData;            % update data
        Re_Plot_Data_Now;      % update plot
        BannedFrames=[BannedFrames;frames2del];
        sprintf('Banned Frames: %d : %d',frames2del)
    end
%% Get NETWROK
%     function Get_Network(~,~)
%        
%     end

end % END OF THE WORLD