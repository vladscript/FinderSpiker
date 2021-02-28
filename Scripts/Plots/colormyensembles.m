% Create Colormaps according to the number
% of conditions and ensembles in each condition
% Using CBREWER function, and some rules in the code for 7 connditions:
% Set1,Dark2,Accent,Set3,Set2,Pastel2,Pastel1
% 
% Input
%   NGroups:    Cell of 1xN conditions and each elemente the scale of
%               N Ensembles in each Condition
% Output
%   ColorState: RGB Matriz of Colors
function ColorState=colormyensembles(NGroups)
%% Setup
NC=numel(NGroups);           % N Conditions
NE=sum(cell2mat( NGroups )); % N Ensembles Total
ColorState=zeros(NE,3);      % Initialize Black Colors
Ncolor=zeros(NE,1);          
%% PRESEST Set COLORS according to CBREWER
%  SEE cbrewer() to see or change colors
fprintf('\n>> Building Colors... ')
% Condition 1 
i=1;
Ncolor(i)=9;
CM=cbrewer('qual','Set1',Ncolor(i));
% First Condition 3 Colors (see Pérez-Ortega et al, 2016)
redensemble=[0.8549,0.1,0.11];
greenensemble=[0.44,0.77,0];
blueensemble=[0.17,0.39,0.99];
CM(1,:)=redensemble;
CM(2,:)=greenensemble;
CM(3,:)=blueensemble;

ConditionColors{i}=CM(setdiff(1:Ncolor(i),4),:); % Without Purples-Kind
Ncolor(i)=8;

% Condition 2
i=2;
Ncolor(i)=8;
CM=cbrewer('qual','Dark2',Ncolor(i));
ConditionColors{i}=CM(end:-1:1,:); % OK colors
% Condition 3
i=3;
Ncolor(i)=8;
CM=cbrewer('qual','Accent',Ncolor(i));
ConditionColors{i}=CM(setdiff(1:Ncolor(i),6),:); % Without Purples-Kind
Ncolor(i)=7;
% Condition 4
i=4;
Ncolor(i)=12;
CM=cbrewer('qual','Set3',Ncolor(i));
ConditionColors{i}=CM(setdiff(1:Ncolor(i),10),:); % Without Purples-Kind
Ncolor(i)=11;
% Condition 5
i=5;
Ncolor(i)=12;
ConditionColors{i}=cbrewer('qual','Set3',Ncolor(i));
% Condition 6
i=6;
Ncolor(i)=8;
ConditionColors{i}=cbrewer('qual','Pastel2',Ncolor(i));
% Condition 7
i=7;
Ncolor(i)=9;
ConditionColors{i}=cbrewer('qual','Pastel1',Ncolor(i));

% Extra Conditoins: Unlikely
if NC>7
    Nensxtr=0;
    for i=8:NC
        Nensxtr=Nensxtr+NGroups{i};
        Ncolor(i)=NGroups{i};
    end
    XtraColors=cbrewer('qual','Paired',Nensxtr,'lineal');
    Nensxtr=1;
    for i=8:NC
        ConditionColors{i}=XtraColors(Nensxtr:Nensxtr+NGroups{i}-1,:);
        Nensxtr=NGroups{i}+1;
    end
end

hubbensemble=[ 0.70588,0,0.78431]; % DEEP PURPLE

fprintf('Done\n')
%% MAIN LOOP: assigning colors:
fprintf('>>Colormap:')
ccounter=1;
for c=1:NC
    if NGroups{c}<Ncolor(c)
        CM=ConditionColors{c}(1:NGroups{c},:);
    else
        CM=cbrewer('qual','Set1',NGroups{c},'lineal');
        CM=CM(end:-1:1,:); % Upside Down
        fprintf('+')
    end
    ColorState(ccounter:ccounter+NGroups{c}-1,:)=CM;
    ccounter=ccounter+NGroups{c};
    fprintf('>#')
end
ColorState=[ColorState;hubbensemble];
fprintf(' Ready.\n')