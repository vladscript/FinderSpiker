% Script that gets the Network from Whole Raster
% Input:
%   R_Condition: Cell of Rasters
%   MetaDataColocaliation: List of [+] & [-] Cells
% Outputs:
% Colored Raster
% CS files to get the Gephi Network
% Update .mat File
%%***************************************************************
% Initial Setup 
Ncells=size(R_Condition{1},1);
StatesofNeurons=ones(Ncells,1);
ColorNeuronCell=zeros(Ncells,3);
C=size(R_Condition,2);
%% SELECT COLOR FOR (+) cells
figureCM=figure;
figureCM.Name='HSV colormap';
figureCM.Position=[612 515 560 118];
imagesc([1:10]);
CM=hsv(10);
figureCM.Colormap=CM;
ColorIndx= inputdlg('Set Color Index [1:10],0 is black, else is pink: ',...
         'Select color for + Cells', [1 70]);
waitfor(ColorIndx);
delete(figureCM);
CIndx=str2num(ColorIndx{1});
if ismember(CIndx,[1:10])
    PositiveCellsColor=255*CM(CIndx,:);
elseif CIndx==0;
    PositiveCellsColor=[0,0,0];
    disp('Color set as deep BLACK')
else
    PositiveCellsColor=255*[1,0.6,0.78];
    disp('Color set as deep PINK')
end
%% Check if There Are Merged Neurons
if exist('MetaDataColocaliation','var')
    aremerged=true;   % Are there already colocated cells
    PositiveCells=MetaDataColocaliation.PositiveCells;
    NegativeCells=MetaDataColocaliation.NegativeCells;
    % Set Color:
    ColorNeuronCell(PositiveCells,:)=repmat(PositiveCellsColor,numel(PositiveCells),1);
else
    aremerged=false;  % Are there already colocated cells
    ColorNeuronCell=repmat(PositiveCellsColor,Ncells,1);
end
%% MAKE char Vector of Colors
ColorString=cell(Ncells,1);
for n=1:Ncells
    for j=1:3
        if j<3
            ColorString{n}=[ColorString{n},num2str(ColorNeuronCell(n,j)),','];
        else
            ColorString{n}=[ColorString{n},num2str(ColorNeuronCell(n,j))];
        end
    end
end
%% Make Networks ****************************************************
for c=1:C
    Raster=R_Condition{c};    
    % TOTAL NETWORK #######################################################
    AdjacencyMatrix=GetAdjacencyMatrix(Raster);
    MaxSynLinks=max(AdjacencyMatrix(:));
    AdjacencyMatrix=AdjacencyMatrix./MaxSynLinks; % NORMALIZED
    % Get Source, Target and Weigth for eache link NO SORTED
    SOURCE=[];
    TARGET=[];
    WEIGHT=[];
    for i=1:Ncells-1
        for j=i+1:Ncells
            if AdjacencyMatrix(i,j)>0
                SOURCE=[SOURCE;i];
                TARGET=[TARGET;j];
                WEIGHT=[WEIGHT;AdjacencyMatrix(i,j)];
            end
        end
    end
    if isempty(SOURCE)
        SOURCE=0;
        TARGET=0;
        WEIGHT=0;
        disp('>>NO NETWORK FOUND')
    end
    % FEATURES ***********************************************
    SynStrengthStats=[mean(WEIGHT),mode(WEIGHT),median(WEIGHT),var(WEIGHT),skewness(WEIGHT),kurtosis(WEIGHT)];
    IndexMaxW=find(WEIGHT==1);
    MaxCoupledPair=[SOURCE(IndexMaxW),TARGET(IndexMaxW)];
    % SAVE CSVs ***********************************************
    FileNameExp = inputdlg(['Save ',Names_Conditions{c},' of ',Experiment,' TOTAL Network as: '],...
             'CSV files for Gephi Visualization', [1 100]);
    if ~isempty(FileNameExp)
        FileNameExp=FileNameExp{:};
            % HEADER:                   State
            % ID: index on Neuron       {OK}
            % Label: Color/States       {OK}
            % Latitude: Y COordinate    {OK}
            % Longitude: X Coordinate   {OK}
            % FileNameExp='Network';
        % NodesTable *******************************************
        CarpetName='\NetWorks-CSV';
        DirSave=pwd;
        poslash=find(DirSave=='\');
        if ~isdir([DirSave(1:poslash(end)-1),CarpetName])
            disp('Folder [NetWorks-CSV] created')
            mkdir([DirSave(1:poslash(end)-1),CarpetName]);
        end
        SaveasName=[CarpetName,'\',Experiment,'_',FileNameExp,'_NODES_TOTAL.csv'];
        HeadersNodes={'ID','Label','Latitude','Longitude','Color'};
        Tnodes=table((1:Ncells)',StatesofNeurons,XY_selected(:,2),XY_selected(:,1),ColorString,...
                'VariableNames',HeadersNodes);
        writetable(Tnodes,[DirSave(1:poslash(end)-1),SaveasName],'Delimiter',',','QuoteStrings',true);

        % Links Table -> K-DEGREE *******************************************
        % Be the linkt between the neuron i and j
        %                                     State
        % Source: Neuron i                      {OK}
        % Target: Neuron j                      {OK}
        % Weight: K-links between them          {OK}
        % Type: Undirected{FIX}                 {OK}
        if numel(SOURCE)>1
            TypeNetwork=repmat('Undirected',numel(SOURCE),1);
        else
            TypeNetwork=mat2cell('Undirected',1);
        end

        SaveasName=[CarpetName,'\',Experiment,'_',FileNameExp,'_EDGES_TOTAL.csv'];
        HeadersLinks={'Source','Target','Weight','Type'};
        % Make Table (without dancing xD)
        Tlinks=table(SOURCE,TARGET,WEIGHT,TypeNetwork,...
            'VariableNames',HeadersLinks);
        writetable(Tlinks,[DirSave(1:poslash(end)),SaveasName],'Delimiter',',','QuoteStrings',true);
        disp('      ***********************');
        disp('      * Saved Total Network *');
        disp('      ***********************');
    else
        disp('[> NETWORK NOT SAVED <]')
    end
    
    
end
if aremerged
    % Plot Raster with separated [+] Cells
    IndexPositiveCells=zeros(Ncells,1);
    IndexPositiveCells(PositiveCells)=1;
    % Plot_Raster_Ensembles(RASTER_Selected_Clean,fs,10,1:Ncells,IndexPositiveCells);
    Plot_Raster_Ensembles(RASTER_Selected_Clean,fs,10,[PositiveCells;NegativeCells],IndexPositiveCells,PositiveCellsColor/255);
    ActualFig=gcf;
    ActualFig.Name=['Highlighted ',MetaDataColocaliation.Cells{1},' Cells'];
end
    