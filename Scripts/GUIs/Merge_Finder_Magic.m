%% Script to Run Colocalization function ----------------------------------
% run afte selection:
% Decalare New Variables
Experiment=Experiment(Experiment~='\');     % NAMES PATCH
global XY_merged;
% global dyename;
% global ColocateIndx;
global MetaDataColocaliation;
[XY_merged,MetaDataColocaliation]=get_merged_coordinates_delux(Experiment,dyename,XY_selected,r);
waitfor(gcf);
%% Update .mat File:
Update_Mergeing;