%% Function to create Venn diagrams for 2 Conditions
% For change in cell activity
% For active or inactive in conditoins
% Input
%   R_Condition Cell of Raster Matrices
%   Names of Conditions: Experimental Conditions Cell String
%       UI asks for pair of conditions : CTRL vs EXP
% Output
%  Set of Nuerons according to the change in its RoA
%  ActiveCells per Pair of Selected Condtions
%  TwoNames Selected COnditions
%  RoA: Vector with RoA for CTRL and EXP
function [ActivityGroups, ActiveCells, TwoNames, RoA]= vennRASTER(R_Condition,Names_Conditions)
%% Init
DeltaTh=0.10; % Threshold to define active, inactive and unchanged
%% Select Pair of conditions
OrdinalString={'1st','2nd'};
for c=1:2
    [index_var{c},index_CHECK{c}] = listdlg('PromptString',['Select ',OrdinalString{c},' Condition'],...
                    'SelectionMode','single',...
                    'Name',['>> Condition: [',num2str(c),'/2] '],...
                    'ListSize',[300 200],...
                    'ListString',Names_Conditions);            
end
%% Get Activity GROUPS
R1=R_Condition{index_var{1}};
R2=R_Condition{index_var{2}};
TwoNames{1}=Names_Conditions{index_var{1}};
TwoNames{2}=Names_Conditions{index_var{2}};
DeltaRoA=getRoA(R1)-getRoA(R2);
ThRoA=max(abs(DeltaRoA))*DeltaTh;
% [px,x]=ksdensity(DeltaRoA,linspace(min(DeltaRoA),max(DeltaRoA),100),'function','pdf');
% % PDF RoA_1-RoA_2
% plot(x,px)
% axis tight; grid on;
% Determind the following groups according to Duhne et al 2020

% Inhibited:    RoA_R1>%Th && RoA_R2=0
ActivityGroups.Ninh= intersect( find(getRoA(R2)==0),find(DeltaRoA>ThRoA) ) ;
% Depressed     RoA_R1>RoA_R2
ActivityGroups.Ndep = find(DeltaRoA>ThRoA);
% Facilitated   RoA_R1<RoA_R2
ActivityGroups.Nfac = find(DeltaRoA<-ThRoA);
% Activated     RoA_R1=0 && RoA_R2>0
ActivityGroups.Nact = intersect( find(getRoA(R1)==0),find(DeltaRoA<-ThRoA) ) ;
% Unchanged     RoA_R1~RoA_R2<%Th
ActivityGroups.Nunc = find(abs(DeltaRoA)<=ThRoA);

RoA.inh=[getRoA(R1(ActivityGroups.Ninh,:)),getRoA(R2(ActivityGroups.Ninh,:))];
RoA.dep=[getRoA(R1(ActivityGroups.Ndep,:)),getRoA(R2(ActivityGroups.Ndep,:))];
RoA.fac=[getRoA(R1(ActivityGroups.Nfac,:)),getRoA(R2(ActivityGroups.Nfac,:))];
RoA.act=[getRoA(R1(ActivityGroups.Nact,:)),getRoA(R2(ActivityGroups.Nact,:))];
RoA.unc=[getRoA(R1(ActivityGroups.Nunc,:)),getRoA(R2(ActivityGroups.Nunc,:))];

%% Get ACTIVE CELLS GROUPS
CellsIndx=1:size(R1,1);
ActiveCells.CellsIndx=CellsIndx;
ActiveCells.AC_1=CellsIndx(sum(R1,2)>0);
ActiveCells.AC_2=CellsIndx(sum(R2,2)>0);
disp('*')

%% Nested Function
function RoA=getRoA(R)
    F=size(R,2);
    RoA=sum(R,2)/F;
end

end
