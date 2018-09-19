%% Script to Run Colocalization function
% run afte selection:
% Decalare New Variables
global XY_merged;
% global ColocateIndx;
global MetaDataColocaliation;
[XY_merged,MetaDataColocaliation]=get_merged_coordinates(Experiment,XY_selected,r);
%% Plot Stuff
[Ntotal,~]=size(XY_selected);
[Nmerged,~]=size(XY_merged);
% Get Indexes of Selected Ones that CoLocolize:
[~,~,MergedIndx]=intersect(XY_merged,XY_selected,'rows');
NotMergedIndx=setdiff(1:Ntotal,MergedIndx);
[~,NC]=size(R_Condition);
% Get Raster from Only Colocated Cells:
for n=1:NC
    R_merged{n}=R_Condition{n}(MergedIndx,:);
    R_nomerged{n}=R_Condition{n}(NotMergedIndx,:);
end
% The POSITIVE ONES ************************************
[Merged_IndexSorted,~,RASTER_MERGED]=SortNeuronsCondition(R_merged);
Plot_Raster_V(RASTER_MERGED(Merged_IndexSorted,:),fs);                           % Clean Whole Raster
set(gcf,'Name',['ID: ',Experiment(2:end),' Merged ',MetaDataColocaliation.Cells{1},' Cells with ',MetaDataColocaliation.Dye{1}],'NumberTitle','off')
Label_Condition_Raster(Names_Conditions,R_merged,fs);   % Labels
% The NEGATIVE ONE(S) **********************************
[Merged_IndexSorted,~,RASTER_NOMERGED]=SortNeuronsCondition(R_nomerged);
Plot_Raster_V(RASTER_NOMERGED(Merged_IndexSorted,:),fs);                           % Clean Whole Raster
set(gcf,'Name',['ID: ',Experiment(2:end),' NO merged '],'NumberTitle','off')
Label_Condition_Raster(Names_Conditions,R_merged,fs);   % Labels
% SELECTION:
[RASTER_Selected_Merged,~,R_MergedOnes,Onsets,New_Index]= Select_Raster_for_NN(fs,R_merged,XY_merged,Names_Conditions,Experiment,0);