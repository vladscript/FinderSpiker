%% Function To save CSV files from Network to Gephi Format
% Functional Network of Neural Ensambles
% Only Significant Activity (!)
% RUN AFTER: Ensemble_Sorting.m

% A�ADIR NODOS DE REFERENCIA ESPACIAL: ESQUINAS DE NEGRO
% A�ADIR NODOS SIN CONEXION PARA AMPLIAR LA DIFERENCIA ENTRE GRUPOS
% Optimizar c�digo 

% Input
%   Raster (frames x cells)     Matrix of Neural Activity
%   XY_selected                 Selected Coordinates (re-sorted)
%   NeuronState {}              Cell of Neurons in each Ensemble
%   Cluster_Indexing            Indexes of Sorted Clustering Ensembles
%   ColorState  *               Color Map for the Ensembles
%   Experiment                  Experiment ID
%   Dialogue Input:             File Name to save

% Structure Ouput
%   fNet.SharedNeuron:          N neurons shared by all ensembles
%   fNet.AM:                    Adjacency Matrix
%   fNet.SynStrengthStats       Weights Link Statistics
%   fNet.MaxSynLinks            Maximum of Links Between Neurons
%   fNet.MaxCoupledPair         Pair of Max Connected Neurons (visualization purpose)
% CSV files 4 Gephi Software:
%   EXPID_Condition_Nodes.csv  ID/Label/Coordinates/State's Colors
%   EXPID_Condition_Links.csv Source/Target/Weight
function fNet=Get_Gephi_Network(Raster,XY_selected,Neurons_State,Cluster_Indexing,ColorState,labels_frames,Experiment)
%% Setup
[~,NStates]=size(Neurons_State); % Total States
ActiveNeurons=[];
for i=1:NStates
    ActiveNeurons=[ActiveNeurons,Neurons_State{i}];
end
% Only Neurons Active in Ensembles and their Coordinates
ActiveNeurons=unique(ActiveNeurons);
% XY_ensambles=XY_selected(ActiveNeurons,:);
XY_ensambles=XY_selected;
% Raster to build functional Network
Raster_ensambles=Raster(:,ActiveNeurons);
[~,C]=size(Raster_ensambles); % N Active Cells 
% Number of Neurons in Ensembles
Ncells=C;
% Number of States
NG=NStates;

%% Ensemble whom Neurons Belong to ********************************
EnsemblesList=unique(labels_frames);    % Ensembles
ColorNeuron=zeros(Ncells,3);            % RGB colors
StatesofNeurons=cell(Ncells,1);
% NSeq=ceil(Ncells/(NG+1));
CountersStates=zeros(1,max(EnsemblesList));
SharedNeruons=0;
% The first Ncells in Cluster_Indexing are the Neuron in Ensembles:
for i=1:Ncells
    StatesN=[];
    % Find for each neuron which ensamble belong to
    for j=1:NG
%         if ismember(ActiveNeurons(i),Neurons_State{j})
        if ismember(Cluster_Indexing(i),Neurons_State{j})
            StatesN=[StatesN,EnsemblesList(j)];
        end
    end
    if isempty(StatesN)                         % No ensemble Neurons
        StatesofNeurons{i}=num2str(0);
        disp('Non-Ensamble Neuron')
    else
        stateword=[];
        for w=1:numel(StatesN)
            if w==numel(StatesN)
                stateword=[stateword,num2str(StatesN(w))];
            else
                stateword=[stateword,[num2str(StatesN(w)),',']];
            end
        end
        StatesofNeurons{i}=stateword;
        % StatesN=num2str( StatesN );
        % StatesofNeurons{i}=StatesN(StatesN~=' ');
    end
    [~,NSin]=size(StatesofNeurons{i}); % N-States neuron is IN
    % Neuron Colors******************************
    %                                           Single-Ensemble Neurons
    if NSin==1              
        Current_State=str2double( StatesofNeurons{i} );
        if Current_State>0
            ColorNeuron(Cluster_Indexing(i),:)=ColorState( Current_State,:);
        else 
            ColorNeuron(Cluster_Indexing(i),:)=[0,0,0]; % Likely Useless
            disp('paint it black')
        end
    %                                           All-Ensembles Neuron    
    elseif NSin==NG 
%     elseif NSin==NG && CountersStates(1)>NSeq
            ColorNeuron(Cluster_Indexing(i),:)=ColorState(end,:); % deep purple
            SharedNeruons=SharedNeruons+1;
            disp('Potential Hub Detected --->')
    %                                           2 or more-ensemble Neurons
    else
        CurretnEnsembles=[];
        % get N chars == ',' +1        
        Nchar=numel(StatesofNeurons{i}==',');
        ComasPos=[0,find(StatesofNeurons{i}==','),Nchar+1];
        for q=1:numel(ComasPos)-1
            prechar=ComasPos(q)+1;
            poschar=ComasPos(q+1)-1;
            CurretnEnsembles=[CurretnEnsembles,str2double( StatesofNeurons{i}(prechar:poschar) )];
        end
        % Choose the ensemble with less neurons:
        [~,minEnsemble]=min(CountersStates(CurretnEnsembles));
        PoorEnsemble=CurretnEnsembles(minEnsemble);
        CountersStates(PoorEnsemble)=CountersStates(PoorEnsemble)+1;
        % [i,PoorEnsemble]
        ColorNeuron(Cluster_Indexing(i),:)=ColorState( PoorEnsemble,:); 
    end
    disp(['Neuron: ',num2str(Cluster_Indexing(i)),' @ Ensemble(s): ',num2str(StatesofNeurons{i})])
end

%% Make Links: ***********************************************
ColorNeuronRGB=255*ColorNeuron;
ColorNeuronCell=cell(Ncells,1);
for i=1:Ncells
    ColorNeuronString=num2str(ColorNeuronRGB(Cluster_Indexing(i),:));
    CharacIndx=setdiff(1:length(ColorNeuronString),find(ColorNeuronString==' '));
    % Save Character
    ColorNeuronStr=[];
    for j=1:length(CharacIndx)-1
        if CharacIndx (j+1)==CharacIndx (j)+1
            ColorNeuronStr=[ColorNeuronStr,ColorNeuronString( CharacIndx(j) )];
        else
            ColorNeuronStr=[ColorNeuronStr,ColorNeuronString( CharacIndx(j) ),','];
        end
    end
    ColorNeuronStr=[ColorNeuronStr,ColorNeuronString( CharacIndx(end) )];
    ColorNeuronCell{i}=ColorNeuronStr;
end
%% Adjacency Matrix:
% AS in Get_Total_Network
% This counts how many frames a Pair of Neurons Fired Together
% Be Cell i-th and Cell j-th
% AdjacencyMatrix(i,j)=Number of join firing;
AdjacencyMatrix=GetAdjacencyMatrix(Raster_ensambles);
MaxSynLinks=max(AdjacencyMatrix(:));
% AdjacencyMatrix=AdjacencyMatrix./MaxSynLinks; % NORMALIZED
% Get Source, Target and Weigth for each link UNSORTED
SOURCE=[];
TARGET=[];
WEIGHT=[];
for i=1:Ncells-1
    for j=i+1:Ncells
        if AdjacencyMatrix(i,j)>0
            SOURCE=[SOURCE;ActiveNeurons(i)];
            TARGET=[TARGET;ActiveNeurons(j)];
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
%% NETWORK OUTPUTS ******************************************************
% N Cells that participate in all Ensembles
fNet.SharedNeruons=SharedNeruons;
% Adjacency (Normalized) Matrix
fNet.AM=AdjacencyMatrix;
% THRESHOLD (?)-> here
% Synaptic Strength: weigths of conecctions statistics
fNet.SynStrengthStats=[mean(WEIGHT),mode(WEIGHT),median(WEIGHT),var(WEIGHT),skewness(WEIGHT),kurtosis(WEIGHT)];
% Maximum Synpatic Strength in Links
fNet.MaxSynLinks=MaxSynLinks;
% Max Coupled Neurons
[~,IndexMaxW]=max(WEIGHT);
fNet.MaxCoupledPair=[SOURCE(IndexMaxW),TARGET(IndexMaxW)];


%% Nodes & Links Table *******************************************************
% Input Dialogues
% Prefix To Save As ...
FileNameExp = inputdlg('Save Network as: ',...
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
    SaveasName=[CarpetName,'\',Experiment,'_',FileNameExp,'_NODES.csv'];
    HeadersNodes={'ID','Label','Latitude','Longitude','Color'};
    Tnodes=table(Cluster_Indexing(1:Ncells)',StatesofNeurons,XY_ensambles(Cluster_Indexing(1:Ncells),2),XY_ensambles(Cluster_Indexing(1:Ncells),1),ColorNeuronCell,...
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
        TypeNetwork=repmat('Undirected',length(SOURCE),1);
    else
        TypeNetwork=mat2cell('Undirected',1);
    end
    
    SaveasName=[CarpetName,'\',Experiment,'_',FileNameExp,'_EDGES.csv'];
    HeadersLinks={'Source','Target','Weight','Type'};
    Tlinks=table(SOURCE,TARGET,WEIGHT,TypeNetwork,...
        'VariableNames',HeadersLinks);
    writetable(Tlinks,[DirSave(1:poslash(end)),SaveasName],'Delimiter',',','QuoteStrings',true);
    disp('      *****************');
    disp('      * Saved Network *');
    disp('      *****************');
else
    disp('>>>> NOT SAVED <<<<<<')
end
%% Check Out Figure
% CatStates=categorical( StatesofNeurons );
% figure; 
% scatter(XY_ensambles(:,1),XY_ensambles(:,2),[],ColorNeuronRGB(ActiveNeurons,:),'filled')
% axis([0,512,0,512]); grid on;
% title('Tissue Location of Active Neurons')
% end