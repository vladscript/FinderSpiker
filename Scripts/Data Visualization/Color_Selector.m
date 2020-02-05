% Function that creates a Color Map for the ExperimentalConditons
% Input
%   Starts a GUI to set the colors
% Ouput
%   CM: Colormap
%   ColorIndx: Struture with colors

function [CM,IndxColor]=Color_Selector(NamesCond)
% Colors: #################################################################
Nconditions=numel(NamesCond);
SetColorMap;
CM=cbrewer(KindMap,ColorMapName,Ncolors);
% figure
figureCM=figure;
figureCM.Name=[ColorMapName,' qualitative colormap from  CBREWER'];
figureCM.Position=[612 515 560 118];
imagesc([1:Ncolors]);
figureCM.Colormap=CM;
figureCM.Children.XTick=1:Ncolors;
figureCM.Children.YTick=[];

% Choose Colors
for n=1:Nconditions
    disp('>> Choose Color:') ;
    ColorIndx{n}= inputdlg(['Set Color for ',...
        NamesCond{n},' Condition'],...
         'Select color for ', [1 35]);
    waitfor(ColorIndx{n});
    IndxColor(n)= str2num( cell2mat(  ColorIndx{n} ));
    if ~ismember(IndxColor(n),1:Ncolors)
        IndxColor(n)=n;
        disp('>>ERROR in the index. Assigned Color :')
        disp(n)
    end
end
delete(figureCM);

fprintf('>>To set other Colormap check: \n>>edit SetColorMap %% script & \n>>run cbrewer() %% function\n')