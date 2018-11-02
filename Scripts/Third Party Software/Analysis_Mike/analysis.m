% raster=leerRaster(); % ONLY ACTIVE NEURONS: Size: Cells x Frames
% plotRaster(raster,1)
% similitudes=corrcoef(raster'); 
similitudes=squareform(1-pdist(raster,'jaccard'));
cluster_index=clusterModularity(similitudes,10000);
% plotClusterRaster(raster,cluster_index,1);
% GpoCo=plotGroupCoactivity(raster,cluster_index);
th=testCoactivityGroup(raster,cluster_index,0.85);
% plotCoactivePeaksByGroup(raster,cluster_index,th);