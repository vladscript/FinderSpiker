% Best fit if power law distribution
% Get the parameters of the best fit of power law distribution.
%
% [slope intercept R2 Chi2 bins] = BestPLfit(data, binsMin, binsMax)
%
% Inputs
% data = data to fit power law distribution.
% binsMin = minimum # of bins.
% binsMax = maximum # of bins.
%
% Output
% slope = parameter of law disttribution
% intercept = parameter of law disttribution
% R2 = coefficient of determination
% Chi2 =
%
% ..:: by Jesús E. Pérez-Ortega ::.. Apr-2012

function [slope intercept R2 Chi2 bins] = BestPLfit(data, binsMin, binsMax)

slopes=zeros(binsMax,1);
intercepts=zeros(binsMax,1);
R2s=zeros(binsMax,1);
Chi2s=zeros(binsMax,1);
for D=binsMin:binsMax
    [y x]=hist(data,D);
    [slope, intercept, R2, Chi2] = logfit(x,y,'loglog');
    %disp(['D:' num2str(D) ' slope:' num2str(slope)...
    %    ' intercept:' num2str(intercept) ' R^2:' num2str(R2)])
    slopes(D)=slope;
    intercepts(D)=intercept;
    R2s(D)=R2;
    Chi2s(D)=Chi2;
end
[R2 bins]=max(R2s(binsMin:binsMax));
bins=bins+binsMin-1;
Chi2=Chi2s(bins);
slope=slopes(bins);
intercept=intercepts(bins);

figure(102)
subplot(4,1,1)
plot(slopes)
title('slopes')
subplot(4,1,2)
plot(intercepts)
title('intercepts')
subplot(4,1,3)
plot(R2s)
title('R2s')
subplot(4,1,4)
plot(Chi2s)
title('Chi2s')