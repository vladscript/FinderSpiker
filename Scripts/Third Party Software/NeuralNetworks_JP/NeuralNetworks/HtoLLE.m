% Plot Hierachical clustering in LLE values

% NNDynA V 5.6 compatible

% JP aug-2014 

% Change Data
DataName='CTR300710B';

States=evalin('base', [DataName '_Analysis.Clustering.VectorStateIndex']);
Peaks=evalin('base',[DataName '_Analysis.Peaks.ActiveDataPeaks']);
Join=evalin('base', [DataName '_Analysis.Peaks.Join']);
GroupColors=evalin('base',[DataName '_Analysis.Colors.States']);

% Join Peaks
PeaksJoin=[];
if (size(Peaks,1)==size(Join,1))
    for i=1:max(Join)
        idx=find(Join==i);
        PeaksJoin(i,:)=sum(Peaks(idx,:),1);
        StatesJoin(i)=States(idx(1));
    end
else
    PeaksJoin=Peaks;
    StatesJoin=States;
end
Groups=max(StatesJoin);

% LLE distance
X=PeaksJoin';
N=size(X,2);
X2=sum(X.^2,1);
Distance=repmat(X2,N,1)+repmat(X2',1,N)-2*X'*X;

if N>50
    Max_N=50;
else
    Max_N=N-1;
end

%% Figure Dimensional Reduction
obj=findobj('name',['Optimal Neighbors (' DataName ')']);
if isempty(obj)
    figure('name',['Optimal Neighbors (' DataName ')'],'numbertitle','off','position',[0 0 600 300]);
else
    figure(obj);
end
clf

warning off
DIs=[];
for Neighbors=2:Max_N
    try
        DataRed=LLE_JP(PeaksJoin',Neighbors,3,Distance);
        DistanceRed=squareform(pdist(DataRed','Euclidean'));
        DI = DunnIdx_JP(Groups,DistanceRed,StatesJoin);
        DIs(Neighbors)=DI;
    end
end

plot(DIs); hold on;
[val NeighborsIDEAL]=max(DIs);
plot(NeighborsIDEAL,val,'*r','MarkerSize',20)

DataRed=LLE_JP(PeaksJoin',NeighborsIDEAL,3,Distance);
DistanceRed=squareform(pdist(DataRed','Euclidean'));
title(['Optimal Neighbors: ' num2str(NeighborsIDEAL) ])
warning on

%% Figure Dimensional Reduction
obj=findobj('name',['Dimensional Reduction (' DataName ')']);
if isempty(obj)
    figure('name',['Dimensional Reduction (' DataName ')'],'numbertitle','off','position',[0 0 600 300]);
else
    figure(obj);
end
clf

% Axes positions
Pos1=[0 0 0.5 1];
Pos2=[0.5 0 0.5 1];

% Colors
BGColor=[1 1 1];
ColorLight=[0.5 0 0];
ColorDark=[0.3 0 0];
ColorTextLines=[0 0 0];

% Plot Data Reduced
axes('outerposition',Pos1)
plot3(DataRed(1,:),DataRed(2,:),DataRed(3,:),'.','color',ColorDark,'MarkerSize',20)
xl=get(gca,'xlim'); yl=get(gca,'ylim'); zl=get(gca,'zlim');
delta=max([abs(xl(1)-xl(2)) abs(yl(1)-yl(2)) abs(zl(1)-zl(2))]);
xl(2)=xl(1)+delta; yl(2)=yl(1)+delta; zl(2)=zl(1)+delta;
xlabel('LLE 1');ylabel('LLE 2');zlabel('LLE 3');
set(gca,'xlim',xl,'ylim',yl,'zlim',zl)
set(gca,'XTick',[],'YTick',[],'ZTick',[])
set(gca,'box','on','view',[0 90])
title('LLE reduction')
set(gca,'Tag',['DimRedAxes1' DataName])

% Plot Data
axes('outerposition',Pos2)
for i=1:Groups
    idx=find(StatesJoin==i);
    plot3(DataRed(1,idx),DataRed(2,idx),DataRed(3,idx),'.','color',GroupColors(i,:),'MarkerSize',40)
    hold on
end
xlabel('LLE 1');ylabel('LLE 2');zlabel('LLE 3');


