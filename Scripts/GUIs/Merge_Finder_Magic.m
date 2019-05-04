%% Script to Run Colocalization function ----------------------------------
% run afte selection:
% Decalare New Variables
Experiment=Experiment(Experiment~='\');     % NAMES PATCH
global XY_merged;
% global dyename;
% global ColocateIndx;
global MetaDataColocaliation;
[XY_merged,MetaDataColocaliation]=get_merged_coordinates_delux(Experiment,dyename,XY_selected,r);
waitfor(gcf);
%% Plot Stuff -------------------------------------------------------------
[Ntotal,~]=size(XY_selected);
[Nmerged,~]=size(XY_merged);
% Get Indexes of Selected Ones that CoLocolize:
[~,~,MergedIndx]=intersect(XY_merged,XY_selected,'rows');
NotMergedIndx=setdiff(1:Ntotal,MergedIndx);
MetaDataColocaliation.PositiveCells=MergedIndx;
MetaDataColocaliation.NegativeCells=NotMergedIndx';
[~,NC]=size(R_Condition);

%% CONDITION LOOP
for n=1:NC
    disp(['Analyzing Rasters: ',num2str(n),'/',num2str(NC)])
    % Sampled Global Raster
    % Get Raster from MERGED (and NOT) Cells:
    R_merged{n}=R_Condition{n}(MergedIndx,:);
    R_nomerged{n}=R_Condition{n}(NotMergedIndx,:);
    
end

%% The POSITIVE ONES ******************************************************
[Merged_IndexSorted,~,RASTER_MERGED]=SortNeuronsCondition(R_merged);
Plot_Raster_Ensembles(RASTER_MERGED(Merged_IndexSorted,:),fs);                           % Clean Whole Raster
set(gcf,'Name',['ID: ',Experiment(2:end),' Merged ',MetaDataColocaliation.Cells{1},' Cells with ',MetaDataColocaliation.Dye{1}],'NumberTitle','off')
Label_Condition_Raster(Names_Conditions,R_merged,fs);   % Labels
% The NEGATIVE ONE(S) **********************************
[Merged_IndexSorted,~,RASTER_NOMERGED]=SortNeuronsCondition(R_nomerged);
Plot_Raster_Ensembles(RASTER_NOMERGED(Merged_IndexSorted,:),fs);                           % Clean Whole Raster
set(gcf,'Name',['ID: ',Experiment(2:end),' NO merged '],'NumberTitle','off')
Label_Condition_Raster(Names_Conditions,R_merged,fs);   % Labels
%% UPDATE MAT FILE *******************************************************
DefaultPath='C:\Users\Vladimir\Documents\Doctorado\Software\GetTransitum\Calcium Imaging Signal Processing\FinderSpiker\Processed Data';
if exist(DefaultPath,'dir')==0
    DefaultPath=pwd; % Current Diretory of MATLAB
end
[FileName,PathName] = uigetfile('*.mat',[' Pick the Analysis File ',Experiment],...
    'MultiSelect', 'off',DefaultPath);
dotindex=find(FileName=='.');
if strcmp(FileName(1:dotindex-1),Experiment)
    checkname=0;
    % SAVE DATA
    save([PathName,FileName],'XY_merged',...
        'R_merged','R_nomerged','MetaDataColocaliation','-append');
    disp([Experiment,'   -> UPDATED (Merged Data)'])
elseif FileName==0
    checkname=0;
    disp('*************DISCARDED************')
else
    disp('Not the same Experiment!')
    disp('Try again!')
end