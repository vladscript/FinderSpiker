%% Script to Run Colocalization function ----------------------------------
% run afte selection:
% Decalare New Variables
Experiment=Experiment(Experiment~='\');     % NAMES PATCH
global XY_merged;
% global ColocateIndx;
global MetaDataColocaliation;
[XY_merged,MetaDataColocaliation]=get_merged_coordinates_delux(Experiment,XY_selected,r);
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
%% Create Figure
PDFsFigure=figure;
% MERGED ONES
h1=subplot(2,4,1); % ITI pdf
h2=subplot(2,4,2); % Length pdf
h3=subplot(2,4,3); % CAG pdf
h4=subplot(2,4,4); % RoA pdf
% NO MERGED ONES
h5=subplot(2,4,5); % ITI pdf
h6=subplot(2,4,6); % Length pdf
h7=subplot(2,4,7); % CAG pdf
h8=subplot(2,4,8); % RoA pdf
title(h1,'InterTransient PDF','FontSize',7)
ylabel(h1,'MERGED')
ylabel(h5,'NO MERGED')
title(h2,'Length Transient PDF','FontSize',7)
title(h3,'CAG PDF','FontSize',7)
title(h4,'RoA PDF','FontSize',7)
hold(h1,'on'); hold(h2,'on'); hold(h3,'on'); hold(h4,'on'); 
hold(h5,'on'); hold(h6,'on'); hold(h7,'on'); hold(h8,'on'); 
%% CONDITION LOOP
for n=1:NC
    disp(['Analyzing Rasters: ',num2str(n),'/',num2str(NC)])
    % Sampled Global Raster
    % Get Raster from MERGED (and NOT) Cells:
    R_merged{n}=R_Condition{n}(MergedIndx,:);
    R_nomerged{n}=R_Condition{n}(NotMergedIndx,:);
    % FEATURE EXTRACTION
    % MERGED ONES
    DATApdfMe=gen_feat_table_merged(R_merged{n},fs,'_MERGED',Experiment,Names_Conditions{n});
    DATApdfNoMe=gen_feat_table_merged(R_nomerged{n},fs,'_NO_MERGED',Experiment,Names_Conditions{n});
    disp('... OK')
    % Plot Merged
    plot(h1,DATApdfMe.isibin,DATApdfMe.isip,'LineWidth',2)
    plot(h2,DATApdfMe.LTbin,DATApdfMe.LTp,'LineWidth',2)
    ksdensity(h3,DATApdfMe.cag);
    ksdensity(h4,DATApdfMe.roa);
    % Plot NO Merged
    plot(h5,DATApdfNoMe.isibin,DATApdfNoMe.isip,'LineWidth',2)
    plot(h6,DATApdfNoMe.LTbin,DATApdfNoMe.LTp,'LineWidth',2)
    ksdensity(h7,DATApdfNoMe.cag);
    ksdensity(h8,DATApdfNoMe.roa);
end
% Fixing Plot
axis(h1,'tight'), axis(h2,'tight'); axis(h3,'tight'); axis(h4,'tight');
grid(h1,'on'); grid(h2,'on'); grid(h3,'on'); grid(h4,'on');
axis(h5,'tight'), axis(h6,'tight'); axis(h7,'tight'); axis(h8,'tight');
grid(h5,'on'); grid(h6,'on'); grid(h7,'on'); grid(h8,'on');
legend(h1,Names_Conditions);
set(h1, 'YScale', 'log'); set(h2, 'YScale', 'log');
set(h5, 'YScale', 'log'); set(h6, 'YScale', 'log');
set(h4, 'XScale', 'log'); set(h8, 'XScale', 'log');
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