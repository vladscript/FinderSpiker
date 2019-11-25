% Make Line Sequence a Vector of Values of Dsicrete Onset Times
% Input
%   Onsets          x values
%   yvalues         y values
%   varargin        if any intervals (optional)
% Output
%   xseq
function [xseq,yseq]=makelineseq(Onsets,yvalues,varargin)
% If there are Intervals
if ~isempty(varargin)
    IntervalsX=varargin{1};
    % Add Initial and Final Points
    Onsets=[min(IntervalsX(:,1));Onsets;max(IntervalsX(:,2))];
    yvalues=[yvalues(1);yvalues;yvalues(end)];
end
xseq=Onsets(1):Onsets(end);
yseq=zeros(size(xseq));
% yseq=[];
% MAIN LOOP: where magic happens:
N=numel(Onsets);
for i=1:N-1
    if i<N-1
        tlin=Onsets(i):Onsets(i+1)-1;
    else
        tlin=Onsets(i):Onsets(i+1);
    end
    % PRE SLOPE
    y2pre=yvalues(i+1);
    y1pre=yvalues(i);
    x2pre=Onsets(i+1);
    x1pre=Onsets(i);
    % SLOPES
    pre_m=(y2pre-y1pre)/(x2pre-x1pre);
    % Y- INTERCEPT
    ytest=yvalues(i);
    xtest=Onsets(i);
    pre_b=ytest-pre_m*xtest;
    % LINE
    ypre=pre_m*tlin+pre_b;
    [~,IndxInterval]=intersect( xseq,tlin);
    yseq(IndxInterval)=ypre;
    %yseq=[yseq;ypre'];
end