% Function that get Lines for the given frames
% It Unifies with Lines the points of a given signal x
function xdeline=getdeline(delineframes,x)
%% Setup ***************************************
xdeline=[];
PointA=1;
xlinsegment=[];
% Make sure there are not smaller segments than 3 samples
delineframes=delineframes(delineframes>3);
delineframes=delineframes(delineframes<numel(x)-3);
N=numel(delineframes);
%% MAIN LOOP ***********************************
for n=1:N+1
    % Set Final Point
    if n<=N
        PointB=delineframes(n);
    else
        PointB=numel(x);
    end
    % Line Parameters
    originsort=x(PointA);
    if n==1
        mslope=0;
        originsort=x(PointB);
    elseif n==N+1
        mslope=0;
    else
        mslope=(x(PointB)-x(PointA))/(PointB-PointA);
    end
    % Make Line
    xlinsegment=mslope*([0:(PointB-PointA)])+originsort;
    if n>1
        xdeline=[xdeline,xlinsegment(2:end)];
    else
        xdeline=[xdeline,xlinsegment];
    end 
    PointA=PointB;
end