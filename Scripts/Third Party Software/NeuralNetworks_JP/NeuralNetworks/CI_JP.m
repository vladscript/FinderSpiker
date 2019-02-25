% Connectivity Index of a Peaks Matrix
% Get connectivity index of neuronal group.
% (This index is an idea of Jesús E. Pérez-Ortega)
%
% [idx] = CI_JP(PM)
%
% Inputs
% PM = binary data as matrix PxC (P = #peaks, C = #cells)
% 
% Outputs
% idx = index of connectivity (1 = All connected; 0 = Nothing connected).
%
% ..:: by Jesús E. Pérez-Ortega ::.. Mar-2012

function [idx] = CI_JP(PM)

[C P]=size(PM);
CM=zeros(C);
for j=1:C-1
    for k=j+1:C
        s=sum(PM(j,:)+PM(k,:)>1);
        CM(k,j)=s; % Matrix of conecctivity
    end
end
CM=CM(2:C,1:C-1);
s=sum(CM,1)'+sum(CM,2)-diag(CM);
realSum=sum(s);
% realDesv=std(s)*(C-1);
maxSum=P*(C-1)^2;       % Maximum of sum
% E1=(realSum-realDesv)/maxSum;   % Evaluation 1
E1=realSum/maxSum; %Mar-2012

idx=E1;
% idx=mean(s);