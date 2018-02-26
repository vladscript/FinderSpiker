% Function to Sort Neurons According to the combination of Ensamble
% Participation
% Input
%   labels_frames: enambles for each frame
%   signif_frames: significant frames (ensamble activation)
%   Experiment: Raster
%   NG: Number of ensembles
% Output
%   New_Order:      Sorted Indexes
%   Neurons_State:  Ensembles Sets
function [New_Order,Neurons_State]=OrderClusters(labels_frames,signif_frames,Experiment,NG)
%% Total Number of Cells
[~,NCells]=size(Experiment);

%% Total possible combinations of gropus - - - - - - - - - - - - - - - 
CO={};
for k=1:NG
    CO{k}=nchoosek(1:NG,k); % combinations for elements
end

%% Neurons in each Ensemble
Neurons_State={};
EnsemblesList=unique(labels_frames);
for g=1:NG
    Current_State=EnsemblesList(g);
    Frames_State=signif_frames(labels_frames==Current_State);
    Neurons_State{g}=find(sum(Experiment(Frames_State,:),1)>0);
end
% ActiveNeurons=unique(Neurons_State);
%% Sorting by Ensemble Combinations
aux=1;
Cummulative_Neurons=[];
Indexes_Group_Comb={};
FirstGroups=[];
NNeurons_Size_Group=[];
for g=1:length(CO)
    [Ncomb,~]=size(CO{g});
    Nneurons_Group=0;
    for gg=1:Ncomb
        % Obtain Neurons ACtive for combiantion of states: CO{g}(gg,:)
        Current_Neurons=unique ( cell2mat( Neurons_State(CO{g}(gg,:)) ) ); 
        ElseGroups=setdiff(1:NG,CO{g}(gg,:));
        AuxCell=Neurons_State(ElseGroups(:)); % Recently Modified
        Else_Neurons=unique (cell2mat( AuxCell ) );
        GroupX=setdiff(Current_Neurons,Else_Neurons);
        Indexes_Group_Comb{aux}=setdiff(GroupX,Cummulative_Neurons);
        Cummulative_Neurons=[Cummulative_Neurons,GroupX];
        FirstGroups(aux)=CO{g}(gg); % Only to order according to: 1,2,3,...
        Nneurons_Group=Nneurons_Group+length(Indexes_Group_Comb{aux});
        aux=aux+1;
    end 
    NNeurons_Size_Group(g)=Nneurons_Group;   % - - - - - > OUTPUT
end

%% Final Sorting
New_Order=[];
for g=1:NG
   Group=find(FirstGroups== FirstGroups(g));
   New_Order=[New_Order,cell2mat( Indexes_Group_Comb(Group) )];
end
Passive_Neurons=setdiff(1:NCells,New_Order);
New_Order=[New_Order,Passive_Neurons]; %  - - - - - - - - - - - - - - -> OUTPUT