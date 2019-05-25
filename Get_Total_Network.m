% Script that gets the Network from Whole Raster
% Run after loading MAT file
% Input:
%   R_Condition: Cell of Rasters
%   MetaDataColocaliation: List of [+] & [-] Cells
% Outputs:
% Colored Raster
% CS files to get the Gephi Network
% Update .mat File
%%***************************************************************
% Initial Setup 
Update_Directory;
Experiment=Experiment(Experiment~='\');
Ncells=size(R_Condition{1},1);
StatesofNeurons=ones(Ncells,1);
ColorNeuronCell=zeros(Ncells,3);
C=size(R_Condition,2);
%% SELECT COLOR FOR (+) cells
Ncolors=12;
ColorMapName='Paired'; % See CBREWER help to see details
figureCM=figure;
figureCM.Name=[ColorMapName,' qualitative colormap from  CBREWER'];
figureCM.Position=[612 515 560 118];
imagesc([1:Ncolors]);
% CM=hsv(10);
CM=cbrewer('qual',ColorMapName,Ncolors);
figureCM.Colormap=CM;
figureCM.Children.XTick=1:Ncolors;
figureCM.Children.YTick=[];
ColorIndx= inputdlg('Set Color Index [1:10], 0 is black (-) Cells, else option will be pink: ',...
         'Select color for + Cells', [1 70]);
waitfor(ColorIndx);
delete(figureCM);
CIndx=str2num(ColorIndx{1});
if ismember(CIndx,[1:Ncolors])
    PositiveCellsColor=255*CM(CIndx,:);
elseif CIndx==0;
    PositiveCellsColor=[0,0,0];
    disp('>>Color set as deep DARKNESS, i.e., black')
else
    PositiveCellsColor=255*[1,0.6,0.78];
    disp('>>Color set as deep PINK')
end
%% Check if There Are Merged Neurons
if exist('MetaDataColocaliation','var')
    disp('>Merge Cells ID: Done')
    aremerged=true;   % Are there already colocated cells
    PositiveCells=MetaDataColocaliation.PositiveCells;
    NegativeCells=MetaDataColocaliation.NegativeCells;
    % For the plotting:
    IndexPositiveCells=zeros(Ncells,1);
    IndexPositiveCells(PositiveCells)=1;
    SortingIndex=[PositiveCells;NegativeCells];
    NameCells=MetaDataColocaliation.Cells{1};
    % Set Color:
    ColorNeuronCell(PositiveCells,:)=repmat(PositiveCellsColor,numel(PositiveCells),1);
else
    disp('>ALL Cells ')
    aremerged=false;  % Are there already colocated cells
    ColorNeuronCell=repmat(PositiveCellsColor,Ncells,1);
    % For the Plotting
    IndexPositiveCells=ones(Ncells,1);
    SortingIndex=1:Ncells;
    NameCells='ALL';
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
NameDir='NetWorks-CSV\';
FileDirSave=pwd;
slashes=find(FileDirSave=='\');
FileDirSave=FileDirSave(1:slashes(end));
step=-0.05;
SetColorMap; % SAME AS AIMs colors
% ColorsMap=cbrewer('qual','Set1',C*2);
ColorsMap=cbrewer(KindMap,ColorMapName,Ncolors);
legobj=[];
for c=1:C
    Raster=R_Condition{c};
    CAG=sum(Raster);
    AUCf=autocorr(CAG);
    AUCcag=AUCf(2);
    Frames=size(Raster,2);
    % TOTAL NETWORK #######################################################
    % As in Get_Gephi_Network***
    AdjacencyMatrix=GetAdjacencyMatrix(Raster);
    MaxSynLinks=max(AdjacencyMatrix(:));
    fprintf('Maximum Time Neurons Connected: %3.2f seconds\n',Frames*MaxSynLinks/fs);
    fprintf('%3.2f %% percentage of the Experiment\n',100*MaxSynLinks);

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
    % NETWORK FEATURES ***********************************************
    SynStrengthStats=[mean(WEIGHT),mode(WEIGHT),median(WEIGHT),var(WEIGHT),skewness(WEIGHT),kurtosis(WEIGHT)];
    [~,IndexMaxW]=max(WEIGHT);
    MaxCoupledPair=[SOURCE(IndexMaxW),TARGET(IndexMaxW)];
    fprintf('Most Linked Neurons: %i <-> %i\n',MaxCoupledPair(1),MaxCoupledPair(2))
    % SAVE CSVs for GEPHI  ***********************************************
    FileNameExp = inputdlg(['Save ',Names_Conditions{c},'Network of ',Experiment,' TOTAL Network as: '],...
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
            disp(['Folder [',CarpetName,'] created'])
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
        disp('[> NETWORK WAS NOT SAVED <]')
    end
    % Save Fetures of [+],[-] OR ALL
    if aremerged
        % Separate [+] & [-] and ALL
        HeadersPopulation={'AUC_pos','AUC_neg',...
            'mean_pp','mode_pp','median_pp','var_pp','skew_pp','kurt_pp',...
            'mean_nn','mode_nn','median_nn','var_nn','skew_nn','kurt_nn',...
            'mean_pn','mode_pn','median_pn','var_pn','skew_pn','kurt_pn'};
        % Positive Cells ***********************************************
        R_pos=R_merged{c};
        CAGpos=sum(R_pos);
        AUCf=autocorr(CAGpos,1);
        AUCcagpos=AUCf(2);
        % AMpos=AdjacencyMatrix(sort(PositiveCells),sort(PositiveCells));
        SOURCEpos=SOURCE(ismember(SOURCE,PositiveCells));
        TARGETpos=TARGET(ismember(SOURCE,PositiveCells));
        WEIGHTpos=WEIGHT(ismember(SOURCE,PositiveCells));
        % Weight by Population Mixtures:
        WPoswPos=WEIGHTpos(ismember(TARGETpos,PositiveCells));
        WPoswNeg=WEIGHTpos(ismember(TARGETpos,NegativeCells));
        % Negative Cells ***********************************************
        R_neg=R_nomerged{c};
        CAGneg=sum(R_neg);
        AUCf=autocorr(CAGneg,1);
        AUCcagneg=AUCf(2);
        % AMneg=AdjacencyMatrix(sort(NegativeCells),sort(NegativeCells));
        SOURCEneg=SOURCE(ismember(SOURCE,NegativeCells));
        TARGETneg=TARGET(ismember(SOURCE,NegativeCells));
        WEIGHTneg=WEIGHT(ismember(SOURCE,NegativeCells));
        % Weight by Population Mixtures:
        WNegwNeg=WEIGHTneg(ismember(TARGETneg,NegativeCells));
        WNegwPos=WEIGHTneg(ismember(TARGETneg,PositiveCells));
        % Merge
        Walte=[WNegwPos;WPoswNeg];
        % Save Stuff
        % Positive Weight Links
        SynStatsPos=[mean(WEIGHTpos),mode(WEIGHTpos),median(WEIGHTpos),var(WEIGHTpos),skewness(WEIGHTpos),kurtosis(WEIGHTpos)];
        % Negative Weight Links
        SynStatsNeg=[mean(WEIGHTneg),mode(WEIGHTneg),median(WEIGHTneg),var(WEIGHTpos),skewness(WEIGHTpos),kurtosis(WEIGHTpos)];
        % Only Positive Weight Links ##################################
        SynPosPos=[mean(WPoswPos),mode(WPoswPos),median(WPoswPos),var(WPoswPos),skewness(WPoswPos),kurtosis(WPoswPos)];
        % Only Negative Weight Links ##################################
        SynNegNeg=[mean(WNegwNeg),mode(WNegwNeg),median(WNegwNeg),var(WNegwNeg),skewness(WNegwNeg),kurtosis(WNegwNeg)];
        % Only Alternate Weight Links #################################
        SynAlt=[mean(Walte),mode(Walte),median(Walte),var(Walte),skewness(Walte),kurtosis(Walte)];
        TfeatsPopulation=table(AUCcagpos,AUCcagneg,...
            SynPosPos(1),SynPosPos(2),SynPosPos(3),SynPosPos(4),SynPosPos(5),SynPosPos(6),...
            SynNegNeg(1),SynNegNeg(2),SynNegNeg(3),SynNegNeg(4),SynNegNeg(5),SynNegNeg(6),...
            SynAlt(1),SynAlt(2),SynAlt(3),SynAlt(4),SynAlt(5),SynAlt(6),...
            'VariableNames',HeadersPopulation);
        % MAKE CSVs
        if ~isdir([FileDirSave,NameDir])
            mkdir([FileDirSave,NameDir]);
            disp('>>Directory >Ensemble Features< created')
        end
        disp('>>Saving...') 
        writetable(TfeatsPopulation,[FileDirSave,NameDir,Experiment,'_',...
            Names_Conditions{c},'_byPOP_fNET','.csv'],...
            'Delimiter',',','QuoteStrings',true);
        disp(['>>Saved at /Ensemble Features: ',Experiment,'-',Names_Conditions{c}])
    end
    % Save ALL cells Features
    %Headers
    HeadersTOTAL={'AUC_tot',...
            'mean_tot','mode_tot','median_tot','var_tot','skew_tot','kurt_tot'};
    %table
    TfeatsALL=table(AUCcag,...
            SynStrengthStats(1),SynStrengthStats(2),SynStrengthStats(3),...
            SynStrengthStats(4),SynStrengthStats(5),SynStrengthStats(6),...
            'VariableNames',HeadersTOTAL);
    %Save
    if ~isdir([FileDirSave,NameDir])
        mkdir([FileDirSave,NameDir]);
        disp('>>Directory >Ensemble Features< created')
    end
    disp('>>Saving...') 
    writetable(TfeatsALL,[FileDirSave,NameDir,Experiment,'_',...
        Names_Conditions{c},'_TOTAL_fNET','.csv'],...
        'Delimiter',',','QuoteStrings',true);
    disp(['>>Saved at /Ensemble Features: ',Experiment,'-',Names_Conditions{c}])
    %% Plot Histogram Plots of the Wieghted Links
    step=step+.2;
    hplot{c}=raincloud_plot(100*WEIGHT,'color',ColorsMap(c,:),'box_on',1,'alphaval',0.5,...
     'box_dodge', 1, 'box_dodge_amount',step , 'dot_dodge_amount', step, 'box_col_match',1,...
     'band_width',0.2);
    legobj=[legobj,hplot{c}{1}];
end
legend(legobj,Names_Conditions);
xlabel('%TimeLinked')
title('PDF')
axis tight; grid on;
FigPDF=gcf;
FigPDF.Name='Distribution of Linking Percentage';
%% Plot Raster ####################################################
Plot_Raster_Ensembles(RASTER_Selected_Clean,fs,1,SortingIndex,IndexPositiveCells,PositiveCellsColor/255);
ActualFig=gcf;
ActualFig.Name=['Highlighted ',NameCells,' Cells of: ',Experiment];
Label_Condition_Raster(Names_Conditions,R_Condition,fs);   % Labels

