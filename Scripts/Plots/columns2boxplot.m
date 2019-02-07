% Function to plot different boxpltos from several column vectors
% Input
%   VC1,VC2,VCN: column vectors 
%   Labels=cell(N,1) Structure of NLabels
% Output 
%   Plot of Several Boxplots (different N data)
function columns2boxplot(varargin)
%% Setup
Nc=numel(varargin);
if iscell( varargin(end) )
    Nboxes=Nc-1;
    fprintf('>> Displaying %i boxplots\n',Nboxes)
    disp('with the Following Labels:')
    Labels=varargin(end);
    disp(cellstr(Labels{1}))
else
    Nboxes=Nc;
    fprintf('>> Displaying %i boxplots',Nboxes)
    for n=1:Nc
        Labels{n,1}=num2str(n);
    end
end
% Create Suitable Data Format
ALL_DATA=[];
Label_Ddata=[];
% TO DO THIS: boxplot([ALL_DATA],[Label_Ddata])
for n=1:Nboxes
    ALL_DATA=[ALL_DATA;varargin{n}];
    Label_Ddata=[Label_Ddata;n*ones(numel(varargin{n}),1)];
end
figure; boxplot(ALL_DATA,Label_Ddata,'Labels',Labels);
grid on;
disp('>> Done.')
    
% boxplot([lDOPA;CLOZA],[ones(6,1);2*ones(4,1)])