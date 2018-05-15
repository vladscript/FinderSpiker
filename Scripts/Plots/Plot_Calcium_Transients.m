% Plot Clean Calcium Transients
% After Manual Data Review
% Input
%   SIGNALSclean,
%   Indexes
%   Onsets
%   R_Condition: Durations
% Output
%   Calcium Transients
function CaTransients=Plot_Calcium_Transients(SIGNALSclean,Onsets,Raster_Condition,Indexes)
% Setup
% New_Index_Active=New_Index(TotalActiveNeurons);
% A2A=find(ColocateIndx);
% Indexes=New_Index_Active(A2A);
NCaT=length(Indexes);
[V,C]=size(SIGNALSclean);

%% Main Loop
CaTransients=[];
for n=1:NCaT
    % Read Signal from Videos and Conditions
    CaTran=[];
    i=Indexes(n);
    % [~,iB]=ismember(XY_selected(Indexes(NCaT),:),XY_clean,'rows');
    for c=1:C % Condition Loop
        Xtrans=[];
        n0=Onsets{c}; % discrete time
        [~,Nsamples]=size(Raster_Condition{c}); % Length
        for v=1:V % Videos Loop
            if ~isempty(SIGNALSclean{v,c})
                xtrans=SIGNALSclean{v,c}(i,:);
                % Normalize
                xtrans(xtrans<0)=0;
                if max(xtrans)>0
                    xtrans=xtrans/max(xtrans);
                end
            end
            Xtrans=[Xtrans,xtrans];
        end
        % Cut Signal
        CaTran=[CaTran,Xtrans(n0:n0+Nsamples-1)];
    end
   CaTransients(n,:)=CaTran;
end
