function CAG_TH_MK=get_CAG_threshold(raster)
rasterAct=raster(find(sum(raster,2)),:);
matriz_adyacente=corrcoef(rasterAct');
[Ci]=clusterModularity(matriz_adyacente,1000);
[significativeThModule]=testCoactivityGroup(rasterAct,Ci,0.9);