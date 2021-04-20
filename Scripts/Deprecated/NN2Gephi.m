%% Function To save CSV files from Network to Gephi Format
% Functional Network of Neural Ensambles
% Only Significant Activity (!)
% RUN AFTER: Ensemble_Sorting.m

% AÑADIR NODOS DE REFERENCIA ESPACIAL: ESQUINAS DE NEGRO
% AÑADIR NODOS SIN CONEXION PARA AMPLIAR LA DIFERENCIA ENTRE GRUPOS
% Optimizar código 

% Input
%   Raster (frames x cells)     Matrix of Neural Activity
%   XY_selected                 Selected Coordinates
%   NeuronState {}              Neurons in each Ensemble
%   Cluster_Indexing            Indexes of Sorted Clustering Ensembles
%   ColorState  *               Color Map for the Ensembles
%   Experiment                  Experiment ID
%   Dialogue Input:             File Name to save

% Ouput
%   SharedNeuron:       N neurons shared by all ensembles
%   Nodes.csv  ID/Label/Coordinates/State's Colors
%   Links.csv Source/Target/Weight
function SharedNeruons=NN2Gephi(Raster,XY_selected,Neurons_State,Cluster_Indexing,ColorState,Experiment)
%% Setup
[~,NStates]=size(Neurons_State);
% ActiveNeurons=[];
% for i=1:NStates
%     ActiveNeurons=[ActiveNeurons,Neurons_State{i}];
% end
% % Only Neurons Active in Ensembles
% ActiveNeurons=unique(ActiveNeurons);
% XY=XY_selected(ActiveNeurons,:);
XY=XY_selected;
% % Raster=Experiment_Analysis.Data.Data;
% Raster=Raster(:,ActiveNeurons);
[frames,C]=size(Raster);
% This counts how many frames a couple of Neurons Fire Together
%% Adjacency Matrix:
% Be Cell i-th and Cell j-th
% AdjacencyMatrix(i,j)=Number of join firing;
AdjacencyMatrix=zeros(C); 
for i=1:frames
            for j=1:C-1         % Cell A
                for k=j+1:C     % Cell B
                    x=Raster(i,j);
                    y=Raster(i,k);
                    if x&&y
                        AdjacencyMatrix(j,k)=AdjacencyMatrix(j,k)+1;
                        AdjacencyMatrix(k,j)=AdjacencyMatrix(k,j)+1;
                    end
                end
            end
end
% Number of Neurons
Ncells=C;
% Number of States
NG=NStates;

%% Write State whom Neurons Belong ********************************
ColorNeuron=zeros(Ncells,3); % RGB colors
StatesofNeurons=cell(Ncells,1);
% NSeq=ceil(Ncells/(NG+1));
CountersStates=zeros(1,NG);
SharedNeruons=0;
for i=1:Ncells
    StatesN=[];
    % Find for each neuron which ensamble belong to
    for j=1:NG
        % if ismember(ActiveNeurons(i),Neurons_State{j})
        if ismember(Cluster_Indexing(i),Neurons_State{j})
            StatesN=[StatesN,j];
        end
    end
    if isempty(StatesN)
        StatesofNeurons{i}=num2str(0);
        disp('Non-Ensamble Neuron')
    else
        StatesN=num2str( StatesN );
        StatesofNeurons{i}=StatesN(StatesN~=' ');
    end
    [~,NSin]=size(StatesofNeurons{i}); % N-States neuron is IN
    % Neuron Colors******************************
    % Single-Ensemble Neurons
    if NSin==1              
        Current_State=str2double( StatesofNeurons{i} );
        if Current_State>0
            ColorNeuron(Cluster_Indexing(i),:)=ColorState( Current_State,:);
        else 
            ColorNeuron(Cluster_Indexing(i),:)=[0,0,0]; % Likely Useless
            disp('paint it black')
        end
        disp('1-ensamble Neuron')
    % All-Ensembles Neuron    
    elseif NSin==NG 
%     elseif NSin==NG && CountersStates(1)>NSeq
            ColorNeuron(Cluster_Indexing(i),:)=ColorState(end,:); % deep purple
            SharedNeruons=SharedNeruons+1;
            disp('Potential Hub Detected --->')
    % 2-ensemble Nuerons
    else
       CurretnEnsembles=[];
        for q=1:numel(StatesofNeurons{i})
            CurretnEnsembles=[CurretnEnsembles,str2double( StatesofNeurons{i}(q) )];
        end
        % Choose the ensemble with less neurons:
        [~,minEnsemble]=min(CountersStates(CurretnEnsembles));
        PoorEnsemble=CurretnEnsembles(minEnsemble);
        CountersStates(PoorEnsemble)=CountersStates(PoorEnsemble)+1;
        % [i,PoorEnsemble]
        ColorNeuron(Cluster_Indexing(i),:)=ColorState( PoorEnsemble,:); 
    end
        
end

%% Make Links: ***********************************************
ColorNeuron=255*ColorNeuron;
ColorNeuronCell=cell(Ncells,1);
for i=1:Ncells
    ColorNeuronString=num2str(ColorNeuron(i,:));
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
% Get Source, Target and Weigth for eache link
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
    
    SaveasName=['\NetWorks-CSV\',Experiment(2:end),'_',FileNameExp,'_Nodes.csv'];
        HeadersNodes={'ID','Label','Latitude','Longitude','Color'};
        Tnodes=table([1:Ncells]',StatesofNeurons,XY(:,2),XY(:,1),ColorNeuronCell,...
            'VariableNames',HeadersNodes);
        DirSave=pwd;
        poslash=find(DirSave=='\');
        writetable(Tnodes,[DirSave(1:poslash(end)-1),SaveasName],'Delimiter',',','QuoteStrings',true);

    % Links Table -> K-DEGREE *******************************************
    % Be the linkt between the neuron i and j
    %                                     State
    % Source: Neuron i                      {OK}
    % Target: Neuron j                      {OK}
    % Weight: K-links between them          {OK}
    % Type: Undirected{FIX}                 {OK}
    TypeNetwork=repmat('Undirected',length(SOURCE),1);
    SaveasName=['NetWorks-CSV\',Experiment(2:end),'_',FileNameExp,'_Links.csv'];
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
% scatter(XY(:,1),XY(:,2),[],ColorNeuron,'filled')
% axis([0,512,0,512]); grid on;
% title('Tissue Location of Active Neurons')
% end