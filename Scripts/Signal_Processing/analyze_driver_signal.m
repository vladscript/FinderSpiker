% Function to anlize Driver Function:
% Check if first samples are the biggest ones (mainly due to miss
% detrending
% Input
%   D:      Original Set of Row Vector of Driver Signals
% Output
%   Dfix:   Fixed Set of Row Vector of Driver Signals
function Dfix=analyze_driver_signal(D)
[C,~]=size(D);
Dfix=D;
for c=1:C
    d=D(c,:);
    % Check amplitudes
    [~,framax]=max(d);
    if framax==1 && sum(d)~=0
        
        nextframe=framax;
        while d(nextframe)>0
            d(nextframe)=0;
            nextframe=nextframe+1;
        end
        Dfix(c,:)=d; % update and fix
        disp('Possible Miss-Detrended Signal: Fixed')
    else
        disp('Driver OK')
    end
end