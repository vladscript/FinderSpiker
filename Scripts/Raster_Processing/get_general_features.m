% Function to extract general features of the raster
% Input
%   R:  Matrix of the Raster of dims: Cell x Frames
% Output
%   Descriptive Features:
%       Descriptive.AN:     Active Neurons
%       Descriptive.eAC     % of Active Cells
%       Descriptive.DurAc   Time of Activity (minutes)
%       Descriptive.CAG     Coactivitygram
%   AI_1:   Mean Activity per Neuron:
%   AI_2:   Effective Activity per Raster

function [Descriptive,AI_1,AI_2]=get_general_features(R)
% General Features:
[C,Frames]=size(R);             % Cells x Frames
AN=sum(sum(R,2)>0);             % Active NEURONS
RAN=AN/C;                       % Ratio of Active Neurons           *
CAG=sum(R);                     % CoActiviGram
CummAc=sum(CAG);                % AUC
DurAc=numel(find(CAG))/Frames;  % Active Frames/ Total Frames
RoA=sum(R')/Frames;             % Rate of Activity per Cell
% INDEXES OF ACTIVITY *********************************************
% Active Ratio CAG
if AN==0
    AI_1=0;
else    
    AI_1=CummAc/(Frames*max(CAG));
end
% Effective Activity
if C==0
    AI_2=0;
else
    AI_2=DurAc*RAN;
end
% OUTPUT ##############################################################
Descriptive.AN=AN;
Descriptive.CAG=CAG;
Descriptive.DurAc=DurAc;
Descriptive.RAN=RAN;
Descriptive.RoA=RoA;
% END [][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][][]