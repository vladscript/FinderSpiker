% Funtion to label multiple data from seveeral experiments of accumulated
% cell groups by drug action 
% data: raw data
% index: vector howmanypointsperexp
% Output
% x: data
% y: experimetn label: '1','2',etc
% 
function [x,y]=getcluslabs(data,index)
Ndata=numel(data);
Nexp=numel(index);
y=cell(Ndata,1);
x=data;
a=1; b=0;
for n=1:Nexp
    b=b+index(n);
    for i=a:b
        y{i}=num2str(n);
    end
    a=b+1;
end
    