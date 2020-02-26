% Gets line from time accumlating function of activity
% Input
%   R: Matrix of activity (cells x frames)
%   fs: sampling frequency
% Output
%   p(1)scalar slope of the line p(2): scalar y-intercept of the line
function p=slopethatraster(R,fs)
% Cummulative sum of coactivitygram
CumCAG=cumsum(sum(R));
xtime=[0:numel(CumCAG)-1]/fs; % in seconds
fprintf('\n>>Fitting: ');
p=polyfit(xtime,CumCAG,1);
fprintf('done.\n');
fprintf('>>Equation y = %3.2fx+%3.2f\nEnd\n',p(1),p(2));