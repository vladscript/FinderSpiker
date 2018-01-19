% Function to Select Raster frames per Condition &
% Extract Descriptive Features.
% Run After Raster_Magic_Better.m (Processing Blocks) or
% Load Saved Data in the EPXEIMENT_DATE.mat Files
% It runs until selection is OK (according to the user)
% Input
%   fs:                 Sampling Frequency
%   Raster_Condition:   Separated Raster in Cells   NOT condition SORTED
%   XY:                 Original Coordinates        NOT condition SORTED
%   Names_Conditions:   Names of Conditions
%   Experiment:         Name of the Experiment
% Ouput
%   Raster per Condition  Cell
%   RASTER_Selected_Clean Matrix
function [RASTER_Selected_Clean,XY_selected,R_Condition]= Select_Raster_for_NN(fs,Raster_Condition,XY,Names_Conditions,Experiment)
okbutton='No';
while ~strcmp('Yes',okbutton)
     %% Show Raster per condition
    [~,NC]=size(Raster_Condition);
    Raster_Selecter={};
    for c=1:NC
        R=Raster_Condition{c};
        Plot_Raster_V(R,fs);
        figure_raster=gcf;
        figure_raster.Name=['Raster ',Names_Conditions{c}];
        % Read Selection
        prompt = {'Start @ min:','Finish @ min:'};
        dlg_title = 'Select Raster ';
        num_lines = 1;
        defaultans = {'0','3'};
        answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
        Minute_Start=str2double(answer{1});
        Minute_End=str2double(answer{2});
        Frame_A=Minute_Start*60*fs+1;
        Frame_B=Minute_End*60*fs;
        OnsetA{c,1}=Frame_A;
        R=R(:,Frame_A:Frame_B);
        Raster_Selecter{c}=R;    
        % Show Selected Raster
        close(figure_raster);
    %     Plot_Raster_V(R,fs);
    %     figure_raster=gcf;
    %     figure_raster.Name=['Raster ',Names_Conditions{c}];
        NVaux{c,1}='1'; % Number of "Videos" or Records
    end
    %% Ignore NonActive Coordinates in the Selection and Sort per Condition
    %[New_Index,Raster_Condition_Sel,RASTER_WHOLE]=SortNeuronsCondition(NVaux,Raster_Selecter);
    [New_Index,Raster_Condition_Sel,RASTER_WHOLE]=SortNeuronsCondition(Raster_Selecter);
    % Plot Whole Sorted Raster:
    RASTER_Selected_Clean=RASTER_WHOLE(New_Index,:);
    XY_selected=XY(New_Index,:);
    % Clean Raster and Coordinates
    ActiveNeurons=find(sum(RASTER_Selected_Clean,2)>0);                 % INDEX of Active NEURONS
    RASTER_Selected_Clean=RASTER_Selected_Clean(ActiveNeurons,:);
    XY_selected=XY_selected(ActiveNeurons,:);                           % Clean Coordinates
    Plot_Raster_V(RASTER_Selected_Clean,fs);                            % Clean Whole Raster
    set(gcf,'Name',['ID: ',Experiment(2:end),' selected '],'NumberTitle','off')
    Label_Condition_Raster(Names_Conditions,Raster_Condition_Sel,fs);       % Labels
    % To Save active and Sorted:
    R_Condition={};
    XY_Condition={};
    AN=[];
    CummAc=[];
%     ActiveNeurons=find(sum(RASTER_Selected_Clean,2)>0);     % Active NEURONS
    TotalN=length(ActiveNeurons);
    for c=1:NC
        % Read Original Raster
        R=Raster_Condition_Sel{c};
        % Sort Raster
        R=R(New_Index,:);
        R_Condition{c}=R(ActiveNeurons,:);
        XY_Condition{c}=XY_selected;
        % Descriptive  Features *********************************************
        AN(c,1)=sum(sum(R,2)>0);          % Active NEURONS
        CummAc(c,1)=sum(sum(R));          % Cummulative Activity
        DurAc(c,1)=length(R)/fs/60;       % Duration [min]
        NameTable{c,1}=Experiment(2:end); % Cell Name
        
    end
    disp([AN',CummAc'])
    % save ina table
    HeadersFeatures={'Experiment','TotalNeurons','Condition','Neurons','CummActivity','Onset','Minutes'};
    Trasterfeatures=table(NameTable,TotalN*ones(NC,1),Names_Conditions,AN,CummAc,OnsetA,DurAc,...
            'VariableNames',HeadersFeatures);
    disp(Trasterfeatures);
    okbutton = questdlg('Selection Alright?');
%     if strcmp('Yes',okbutton) % Save Table
%         NameFile=[pwd,'\Results\Rasters_Descriptive_Features.xls'];
%         disp('Selection: done')
%         if exist(NameFile,'file')>0 % If the File Already Exists
%             disp('Writing Table ...')
%             Tprevious=readtable(NameFile);
%             [ColumnsUsed,~]=size(Tprevious);
%             corner_indx=['A',num2str(ColumnsUsed+2)];
%             writetable(Trasterfeatures,NameFile,'Sheet',1,'Range',corner_indx,...
%                 'WriteVariableNames',0);
%         %% For the First Time The Program Save Information
%         else                            % New File
%             disp('Writing NEW Table ...')
%             writetable(Trasterfeatures,NameFile,'Sheet',1);
%         end
%         disp('                       Saved Data ');
%     else
%         disp('Select again ...')
%         
%     end
end
%% Save Rasters, Coordiantes & Features
% RASTER_Selected_Clean:    Matrix
% XY_selected:              matrix
% R_Condition:              cell
% XY_Condition:             cell
