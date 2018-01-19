%% Function to order Neurons by Condition (influented by Plot_Clusterin_from_NN.m)
% Input
%   NumberofVideos:     Number of Videos for Each Condition (vector)
%   Raster:             Raster per viedo per condtion (cell)
% Output
%   NewIndex:           Sorted Neurons by Condition
%   Raster_Conditon:    Raster per condition only
function [New_Index,Raster_Condition,RASTER_WHOLE]=SortNeuronsCondition(RASTER)
% Size of the Cell:
[NV,NC]=size(RASTER);
% RASTER dimensions: Cells x Frames
[NCells,Nframes]=size(RASTER{1});
if NCells>Nframes
    %aux1=NCells;
    NCells=Nframes;
    %Nframes=aux1;
    TrnpseInd=1;
else
    TrnpseInd=0;
end
%% Total possible combinations of conditions - - - - - - - - - - - - - - - 
CO={};
for k=1:NC
    CO{k}=nchoosek(1:NC,k); % combinations for elements
end
%% Active Neurons in each Condition & Reading Rasters
NeuronsCondition={};
RASTER_WHOLE=[];
for i=1:NC
    Raster=[];
    for j=1:NV
        R=RASTER{j,i};  % Read Raster
        if TrnpseInd R=R'; end
        if ~isempty(R)
            Raster=[Raster,R]; % Cummulative Raster
        end
    end
    Raster_Condition{i}=Raster;
    NeuronsCondition{i}=find(sum(Raster,2)>0);
    RASTER_WHOLE=[RASTER_WHOLE,Raster];
end
%% Sorting Magic

aux=1;
Cummulative_Neurons=[];
Indexes_Condition_Comb={};
FirstGroups=[];
NNeurons_Size_Group=[];
for g=1:length(CO)
    [Ncomb,~]=size(CO{g});
    Nneurons_Group=0;
    for gg=1:Ncomb
        % Obtain Neurons ACtive for combiantion of states: CO{g}(gg,:)
        AuxCell=NeuronsCondition(CO{g}(gg,:));
        Current_Neurons=unique ( cell2mat( AuxCell(:) ) ); 
        ElseGroups=setdiff(1:NC,CO{g}(gg,:));
        AuxCell=NeuronsCondition(ElseGroups);
        Else_Neurons=unique (cell2mat( AuxCell(:) ) );
        
        GroupX=setdiff(Current_Neurons,Else_Neurons);
        if ~isempty(GroupX)
            AuxNeurons=setdiff(GroupX,Cummulative_Neurons);
            if ~isempty(AuxNeurons)
                Indexes_Condition_Comb{aux}=AuxNeurons;
            else
                Indexes_Condition_Comb{aux}=[];
            end
        else
            Indexes_Condition_Comb{aux}=[];
        end
        Cummulative_Neurons=unique([Cummulative_Neurons;GroupX]);
        FirstGroups(aux)=CO{g}(gg); % Only to order according to: 1,2,3,...
        Nneurons_Group=Nneurons_Group+length(Indexes_Condition_Comb{aux});
        aux=aux+1;
    end 
    NNeurons_Size_Group(g)=Nneurons_Group;   % - - - - - > OUTPUT
end
% Bar of N neurons in combination of M groups * * * * * * *
% figure; bar(1:NC,NNeurons_Size_Group); grid on; axis tight;
% ylabel('N-Neurons'); xlabel('Belong to M Groups')

New_Index=[];
for g=1:NC
   Group=find(FirstGroups== FirstGroups(g));
   AuxCell=Indexes_Condition_Comb(Group);
   New_Index=[New_Index;cell2mat( AuxCell(:) )];
end
Passive_Neurons=setdiff(1:NCells,New_Index);
New_Index=[New_Index;Passive_Neurons'];