function showmeyourraster(R,Indexes,fs,NSECONDS,NameOutRaster)
fprintf('>>Raster Animation Starting:\n')
v = VideoWriter([NameOutRaster,'.avi']);
% Plot Raster
Plot_Raster_Ensembles(R(Indexes,:),fs,1,[1:numel(Indexes)])
FigRaster=gcf;
%CAGAxis=FigRaster.Children(1);
RasterAxis=FigRaster.Children(2);
% Showtime raster
RasterAxis.XLim=[-NSECONDS/60,0];
F=size(R,2);
DeltaMin=1/(60*fs);
open(v);
for f=1:F
    RasterAxis.XLim=[-NSECONDS/60+f*DeltaMin,0+f*DeltaMin];
    drawnow;
    Fmovie=getframe(gcf);
    writeVideo(v,Fmovie);
    fprintf('Video Animation %3.2f %%\n',100*f/F)
end
close(v); close(gcf);
fprintf('>>Done.\n')