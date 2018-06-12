%% Function that makes Matrix: Cells x Frames
% Assumption:
% Frames>Cells
% Input
%   R: matrix size: DxG
% Output
%  R: maatrix size: 
%  if D>G
function Rout=makecellsxframes(Rin)
[D,G]=size(Rin);
if D>G
    Rout=Rin'; % TRANPOSE
else
    Rout=Rin;
end
