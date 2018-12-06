%%
Plot_Raster_V(RASTER_Selected_Clean)
%%
CAG_TH=4;
R_Dyskinesia=R_Condition{1};
R_Clozapine=R_Condition{2};
R_Dysk_Analysis = raster_cluster(R_Dyskinesia,CAG_TH,3,'hamming');
R_Clz_Analysis = raster_cluster(R_Clozapine,CAG_TH,3,'hamming');
Ensemble_Sorting
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>END')
%%
[Features_Ensemble,Features_Condition]=get_ensembles_features(R_Condition,Ensemble_Threshold,Ensembled_Labels,fs);
save_features_ensembles(Experiment,Names_Conditions,Features_Ensemble,Features_Condition);