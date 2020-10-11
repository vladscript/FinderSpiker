% Gets line from time accumlating function of activity
% 
% CumCAG=p(1)+p(2)*T+W
% Where:
% p: linear fit coeficients
% T: time in seconds [s]
% W: residue
% 
% Input
%   R: Matrix of activity (cells x frames)
%   fs: sampling frequency
% 
% Output
%   p(1)scalar slope of the line p(2): scalar y-intercept of the line
function [p,CumCAG]=slopethatraster(R,fs,varargin)
if isempty(varargin)
    TotalAct=1;
else
    TotalAct=varargin{1};
end
% Cummulative sum of coactivitygram
CumCAG=cumsum(sum(R));
% Normalize to 100% as the max of cumCAG:
PercProportion=CumCAG(end)/TotalAct;
CumCAG=PercProportion*100*CumCAG/CumCAG(end);

% WHY NOT normalize:
% Supose CAG1/10=CAG2
% The same will hava

% Time in Seconds
xtime=[0:numel(CumCAG)-1]/fs; % in seconds
fprintf('\n>>Fitting: ');
p=polyfit(xtime,CumCAG,1);
% %% Check  p(0) make closer to Zero, reducing
% StepPerc=5;
% nSteps=1;
% while abs(p(2))>20
%     % Probably activity raised up rapidly.
%     % Soluction: reduci time interval
%     if nSteps<2 
%         fprintf('\n>>Fixing Fit Interval:')
%     end
%     NsampeOK=find(CumCAG<=100-nSteps*StepPerc);
%     p=polyfit(xtime(NsampeOK),CumCAG(NsampeOK),1);
%     nSteps=nSteps+1;
%     fprintf(' * ')
% end
% fprintf('\n')
fprintf('done.\n');
fprintf('>>Equation y = %3.2fx+%3.2f\nEnd\n',p(1),p(2));