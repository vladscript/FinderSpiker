% This function Retrieves inter events interval & event duration 
% Input
%   r:   vector of ones that indicate the events
% Output
%   IEI: inter events intervals
%   ED:  events duration
function [IEI,ED]=interval_duration_events(r)
% Setup
IEI=[];
ED=[];
if ~isempty(r(r>0))
    dr=diff(r); 
    % +dr -> start Interval
    % -dr -> end Interval
    StartsIntervasl=find(dr>0);
    EndsIntervasl=find(dr<0);
    if and(~isempty(StartsIntervasl),~isempty(EndsIntervasl))
        % If it started at Initial Transient
        if StartsIntervasl(1)>EndsIntervasl(1)
            EndsIntervasl=EndsIntervasl(2:end);
        end
        % If it ended at with out End of Transient
        if StartsIntervasl(end)>EndsIntervasl(end)
            StartsIntervasl=StartsIntervasl(1:end-1);
        end
        if numel(StartsIntervasl)==numel(EndsIntervasl)
            if numel(StartsIntervasl)>1
                IEI=[IEI,StartsIntervasl(2:end)-EndsIntervasl(1:end-1)+1];
            end
            ED=[ED, EndsIntervasl-StartsIntervasl+1];
        else
            disp('WTF?');
        end
    end
else
    disp('>>No events found.')
end