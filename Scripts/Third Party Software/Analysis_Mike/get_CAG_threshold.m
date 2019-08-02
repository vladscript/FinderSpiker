function Ci=get_CAG_threshold(raster,pval)
rasterAct=raster(find(sum(raster,2)),:);
matriz_adyacente=corrcoef(rasterAct');
[Ci]=clusterModularity(matriz_adyacente,1000);
[CAG_TH_MK]=testCoactivityGroup(rasterAct,Ci,pval);