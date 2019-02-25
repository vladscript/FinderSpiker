% Raster with only Active Cells
%
% Cells that are no active are eliminated from raster.
%
% [Xa NoActive Active] = RActive_JP(X, show)
%
% Inputs
% X = binary data as matrix FxC (F = #frames, C = #cells)
% show = show index of no active cells (1=yes; 0=no) 
% 
% Outputs
% Xa = binary data as matrix FxC with only active cells.
% NoActive = indices of no active cells.
% Active = indices of active cells.
%
% ..:: by Jesús E. Pérez-Ortega ::.. Mar-2012

function [Xa NoActive Active] = RActive_JP(X,show)

Cells=sum(X,1);
Active=find(Cells);
NoActive=find(Cells==0);
Xa=X(:,Active);

if show
    disp(['Active Cells: ' num2str(Active)])
    disp(['No Active Cells: ' num2str(NoActive)])
end
