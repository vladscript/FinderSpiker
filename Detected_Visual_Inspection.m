% # Check False Positive 
% Accepted and processed Data
global indxSIGNALSOK;
indxSIGNALSOK=cell(size(isSIGNALS));
indxSIGNALSOK = Calcium_Magic(isSIGNALS);
VisualInspector(1)=true;