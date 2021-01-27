% Make binary matrix from categorical vector
% Input
%   Y:      categorical vector
% Output
%   Ybin:   Binary matrix from Y
% Example
% >>unique(Y)=['Apple';'Pie'];
% >>unique(Ybin,'rows')=[0,1;1,0];
function Ybin=getbinarylabel(Y)
Nobservations=numel(Y);
Labels=unique(Y,'legacy');
Nclasses=numel(Labels);
Ybin=zeros(Nobservations,Nclasses);
for i=1:Nobservations
    Ybin(i,Labels==Y(i))=1;
end
