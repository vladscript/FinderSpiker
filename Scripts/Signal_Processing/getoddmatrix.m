% Funtion to get Odd Matrix from
% Input
% isSIGNALS:        detected signals with Ca++ Transients Automatically
% notSIGNALS:       retected signals Automatically
% Global Input
% indxSIGNALSOK:    true signals with Ca++ Transients
% SIGNALS:          Original Raw Data
% Output
% OddsMatrix: Cell Of Structure Matrix for ith video and ith condition:
%       OddsMatrix{j,i}.TruePositive=TruePositive;          % ++
%       OddsMatrix{j,i}.FalsePositive=FalsePositive;        % -+
%       OddsMatrix{j,i}.FalseNegative=[];                   % -- 
%       OddsMatrix{j,i}.TrueNegative=[];                    % +-
function OddsMatrix=getoddmatrix(isSIGNALS,notSIGNALS)
% Call Global Variable to determine if there're data
global SIGNALS;
global indxSIGNALSOK;
% Setup
[NV,NC]=size(isSIGNALS);
% Initialize Output
OddsMatrix=cell(NV,NC); % Odds Matrix: {++,-+,+-,--}
[Cells,~]=size(SIGNALS{1});
for j=1:NV
    for i=1:NC
        if ~isempty(SIGNALS{j,i})
            indxSIGNALSNOT=setdiff(1:Cells,indxSIGNALSOK{j,i});
            OddsMatrix{j,i}.TruePositive=intersect(indxSIGNALSOK{j,i},isSIGNALS{j,i});
            OddsMatrix{j,i}.FalsePositive=setdiff(indxSIGNALSOK{j,i},isSIGNALS{j,i});
            OddsMatrix{j,i}.FalseNegative=setdiff(indxSIGNALSNOT,notSIGNALS{j,i});
            OddsMatrix{j,i}.TrueNegative=intersect(indxSIGNALSNOT,notSIGNALS{j,i});
        end
    end
end