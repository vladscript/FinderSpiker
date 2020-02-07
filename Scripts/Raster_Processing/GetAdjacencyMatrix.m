% Get Adjacency Matrix from Activity Matrix Normalized by Time (Raster)
% Fire together Wire together
% Assuming having more frames than cells ALWAYS!
% Input
%   R: Raster dim: Cells (in CLSUTER ANALYSIS) x Frames 
% Output
%   A: Adjacency Matrix dim: Cells x Cells where:
%       a_ij: percentage of time that the i-th & j-rh fired together
function AdjacencyMatrix=GetAdjacencyMatrix(Raster)
[C,frames]=size(Raster);
% It is supossed to have more frames than cells ALWAYS!
if frames<C
    Raster=Raster';     % Tranpose Raster
    bC=C;               % backup variable
    C=frames;           % Switch Values
    frames=bC;          % Switch Values
    disp('> Data Transposed');
end
% This counts how many frames a couple of Neurons Fire Together
AdjacencyMatrix=zeros(C); 
ActiveFrames=find(sum(Raster));
for f=1:length(ActiveFrames)% Frames loop
    actualframe=ActiveFrames(f);
    ActiveNeurons=find(Raster(:,actualframe));
    AdjacencyMatrix(ActiveNeurons,ActiveNeurons)=AdjacencyMatrix(ActiveNeurons,ActiveNeurons)+1;
end
% Set Zero to the DIAGONAL:
AdjacencyMatrix=AdjacencyMatrix.*~eye(size(AdjacencyMatrix));
% MaxSynLinks=max(AdjacencyMatrix(:));
% AdjacencyMatrix=AdjacencyMatrix./frames; % NORMALIZED