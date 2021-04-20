% Finder Spiker for Calcium Imaging
% Initial interface to select
%   - Multiple-Condition Experiment
%   - - Detect Spikes from original .AVI videos and
%   - - List of Coordinates (ImPatch) or
%   - - ROI coordinates (ImageJ)
%   - Single-Condition Experiment
%   - - Detect Spikes from Fluorescence Traces (.CSV)
%% Standard Initialization
clc; clear; close all;
Import_FinderSpiker;
%% Show Menu
AnalysisMode=SelectMode();
%% Start Analysis:
switch AnalysisMode
    case'Multiple'
        fprintf('Calcium Imaging:Impatch Analysis\n');
        Multiple_Analysis;
	case'Single'
        fprintf('Calcium Imaging:CSV file Analysis\n')
        Single_Analysis;
    otherwise
        fprintf('\nNo analysis mode selected\nBye!\n')
end
    