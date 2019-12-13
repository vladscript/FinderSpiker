% Funtion to get oval pixels in retangular ROIs
% Input
%   sROI: Cell of ROIs from ImageJ
% Output
%   XYcenter: Coordinates of ROI center
%   ROIpixel: pixels Coordinates in Elliptical ROI
function [XYcenter,ROIpixel]=getpixelsatROI(sROI)
N= numel(sROI);
XYcenter=[];
for n=1:N
    % Get Rectangular Coordinates
    xyrect=sROI{n}.vnRectBounds;
    x1=xyrect(1);
    x2=xyrect(3);
    y1=xyrect(2);
    y2=xyrect(4);
    radiusX=abs(x1-x2)/2;
    radiusY=abs(y1-y2)/2;
    centerX=min([x1,x2])+abs(x1-x2)/2;
    centerY=min([y1,y2])+abs(y1-y2)/2;
    XYcenter=[XYcenter;centerX,centerY];
    % Pixel Mesh
    [columnsInImage rowsInImage] = meshgrid(x1:x2, y1:y2);
    ellipsePixels = (rowsInImage - centerY).^2 ./ radiusY^2 ...
    + (columnsInImage - centerX).^2 ./ radiusX^2 <= 1;
    ROIpixel{n}=[columnsInImage(find(ellipsePixels)),rowsInImage(find(ellipsePixels))];
end
    
