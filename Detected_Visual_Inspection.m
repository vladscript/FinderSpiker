%% Setup
Import_FinderSpiker;
%% # Check False Positive 
% Accepted and processed Data
global indxSIGNALSOK;
indxSIGNALSOK=cell(size(isSIGNALS));

[C,V]=size(SIGNALS);
Ndata=0;
Nvidok=0;
for c=1:C
    for v=1:V
        if ~isempty(SIGNALS{c,v})
            Ndata=Ndata+1;
            if ~isempty(isSIGNALS{c,v})
                Nvidok=Nvidok+1;
            end
        end
    end
end

if Nvidok>0
    indxSIGNALSOK = Calcium_Magic(isSIGNALS);
    VisualInspector(1)=true;
else
    disp('>>NO DETECTED CELLS: all cells were undetected')
end