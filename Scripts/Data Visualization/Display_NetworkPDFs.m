% Requires to load myNetworkDataSet_XXXX.mat FILE
% Directory:
clear; %close all;
Dirpwd=pwd;
slashesindx=find(Dirpwd=='\');
CurrentPathOK=[Dirpwd(1:slashesindx(end))];
[FileName,PathName] = uigetfile('.mat','Load DATASET of GEPHI Features',...
    CurrentPathOK);
fullFileName = fullfile(PathName, FileName);
load(fullFileName);
% Plot P(x) and c(k) from Experiment Networks
[Ne,Nc]=size(GEPHIDATA); % Data are strings
% Choose Colors
[CM,ColorIndx]=Color_Selector(Names_Conditions);

%% Select Parameter
X=GEPHIDATA{1,1};
AllParams=X.Properties.VariableNames;

[index_var,index_CHECK] = listdlg('PromptString','Select Parameter',...
                    'SelectionMode','single',...
                    'Name','Network Parameters',...
                    'ListSize',[300 200],...
                    'ListString',AllParams(5:end));

Ncol=index_var+5-1;
NetParam=AllParams(Ncol);
%% DATA
% EDIT to load data from functions: Display_NetworkPDFs !!!!!!!!!!!!

[TotalTable,DATAnet]=TableNetwork(ExsitFeature,GEPHIDATA,EXPLIST,Names_Conditions,ActualFeature);
[ExsitFeature,EXPLIST]=emptycolnettables(GEPHIDATA,Names_Conditions,ActualFeature);


%% Do Magic Plots

for e=1:Ne
    X=GEPHIDATA{e,:};
    figure('Name',sprintf('CDFs Experiment: %s ',X.EXP_ID(1,:)),...
        'NumberTitle','off','Position',[20 246 1336 420]);
    ECDGaxis{e}=subplot(1,2,1);
    sCDGaxis{e}=subplot(1,2,2);
    figure('Name',sprintf('Experiment: %s ',X.EXP_ID(1,:)),...
        'NumberTitle','off','Position',[20 246 1336 420]);
    MaxK=0;
    MaxParam=0;
    PDF_k=cell(1,Nc);
    PDF_param=cell(1,Nc);
    Nbins_k=zeros(1,Nc);
    Nbins_param=zeros(1,Nc);
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % Use Nodes with Degree>0 in C_i or peviously Degree>0 in C_i-1
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    preActive=[]; 
    for c=1:Nc
        DataNet=GEPHIDATA{e,c};
        Ndata=numel(DataNet.degree);
        var_k=zeros(Ndata,1);
        netfeat=zeros(Ndata,1);
        % Degree ***************************************************
        for n=1:Ndata
            var_k(str2double( DataNet.id{n}),1) = str2double( DataNet.degree{n});
            netfeat(str2double( DataNet.id{n}),1) = str2double( DataNet{:,Ncol}{n} );
        end
        % Current or Previously Active NODES
        preActive=unique([preActive;find(var_k>0)]);
        var_k=var_k(preActive);
        netfeat=netfeat(preActive);
        % Range Values
        if max(var_k)>0
            k_range=1:max(var_k);
        else
            k_range=unique(var_k);
        end
        % Var to make boxplots:
        NetParametersExps{e,c}=netfeat;
        
        MaxK=max([var_k;MaxK]);
        MaxParam=max([netfeat;MaxParam]);
        
        % CDFs: empmirical & smoothed
        [ecdFun,epnet] = ecdf(netfeat);
        [scdFun,spnet]=ksdensity(netfeat,linspace(min(netfeat),max(netfeat),numel(ecdFun)),'function','cdf');
        plot(ECDGaxis{e},epnet,ecdFun,'LineWidth',2,'Color',CM(ColorIndx(c),:));
        hold(ECDGaxis{e},'on')
        plot(sCDGaxis{e},spnet,scdFun,'LineWidth',2,'Color',CM(ColorIndx(c),:));
        hold(sCDGaxis{e},'on')
% %         % KSDENSITY METHOD
% %         % P(k)
% %         [pk,kbin]=ksdensity(var_k,linspace(min(var_k),max(var_k),10),'function','pdf');
% %         % P(parameter)
% %         [pNet,Netbin]=ksdensity(netfeat,linspace(min(netfeat),max(netfeat),10),'function','pdf');
        % HISTOGRAM BASED
%         [pk,edges] = histcounts(var_k,'BinMethod','fd','Normalization','pdf');
        %'auto' | 'scott' | 'fd' | 'integers' | 'sturges' | 'sqrt'
        % Normalize Degree
        k_norm=k_range/k_range(end);
        
        % Parameter as function of the degree
        scaterx=[];
        scatery=[];
        for k=1:k_range(end)
            scatery=[scatery;netfeat(var_k==k)];
            scaterx=[scaterx;k*ones(size(netfeat(var_k==k)))];
        end
        %
        % PLOTTING
        
        % Degree Histogram
        subplot(1,3,1)
        % plot(kbin,pk,'o','Color',CM(ColorIndx(c),:)); hold on;
        % bar(kbin,pk,'FaceColor',CM(ColorIndx(c),:),'EdgeColor','none','BarWidth',0.4)
        PDF_k{c}=histogram(var_k,'BinMethod','integer','Normalization',...
            'probability','FaceAlpha',0.4,'EdgeColor','none','FaceColor',CM(ColorIndx(c),:));
        hold on;
        xlabel('k')
        ylabel('p(k)')
        axis([0,MaxK,0,1]); grid on;
        title('Degree PDF')
        Nbins_k(c)=PDF_k{c}.NumBins;
        % Parameter histogram
        subplot(1,3,2)
        % plot(Netbin,pNet,'o','Color',CM(ColorIndx(c),:)); hold on;
        % bar(Netbin,pNet,'FaceColor',CM(ColorIndx(c),:),'EdgeColor','none','BarWidth',0.4)
        PDF_param{c}=histogram(netfeat,'BinMethod','sturges','Normalization',...
            'probability','FaceAlpha',0.4,'EdgeColor','none','FaceColor',CM(ColorIndx(c),:));
        hold on;
        xlabel(NetParam{1})
        ylabel(sprintf('p(%s)',NetParam{1}))
        axis([0,MaxParam,0,1]); grid on;
        title(sprintf('%s PDF',NetParam{1}))
        Nbins_param(c)=PDF_param{c}.NumBins;
        % Parameter as function of degree
        subplot(1,3,3)
        scatter(scaterx,scatery,[],CM(ColorIndx(c),:),'filled'); hold on;
        xlabel('k')
        ylabel(sprintf('%s',NetParam{1}))
        axis([0,MaxK,0,MaxParam]); grid on;
        title(sprintf('%s (k)',NetParam{1}))
        
%         histogram(var_k,'BinMethod','fd','FaceColor',CM(ColorIndx(c),:));
%         hold on; axis tight; grid on;
    end
    % Fixing Bining
    [~,BestbinK]=max(Nbins_k);
    [~,BestbinParam]=max(Nbins_param);
    ChangeHists_k=setdiff([1:Nc],BestbinK);
    ChangeHists_param=setdiff([1:Nc],BestbinParam);
    for c=1:Nc-1
        PDF_k{ChangeHists_k(c)}.BinEdges=PDF_k{BestbinK}.BinEdges;
        PDF_param{ChangeHists_param(c)}.BinEdges=PDF_param{BestbinParam}.BinEdges;
    end
    for c=1:Nc
        subplot(1,3,1); hold on;
        pk=PDF_k{c}.Values;
        edges=PDF_k{c}.BinEdges;
        plot(edges(1:end-1)+diff(edges)/2,pk,'LineWidth',2,...
            'Color',CM(ColorIndx(c),:),'Marker','o','MarkerSize',4,...
            'MarkerEdgeColor',CM(ColorIndx(c),:),'MarkerFaceColor',CM(ColorIndx(c),:));
        subplot(1,3,2); hold on;
        pParam=PDF_param{c}.Values;
        edgesParam=PDF_param{c}.BinEdges;
        plot(edgesParam(1:end-1)+diff(edgesParam)/2,pParam,'LineWidth',2,...
            'Color',CM(ColorIndx(c),:),'Marker','o','MarkerSize',4,...
            'MarkerEdgeColor',CM(ColorIndx(c),:),'MarkerFaceColor',CM(ColorIndx(c),:));
    end
end
disp('ok')

%% BOXPLOTS
for e=1:Ne
    NetValues=cell(Nc,1);
    for c=1:Nc
        NetValues{c}=NetParametersExps{e,c};
    end
    figure; pairedraincloud('unpaired',CM(ColorIndx,:),0,0,NetValues);
    title(sprintf('%s boxplots of %s',NetParam{1},GEPHIDATA{e,1}.EXP_ID(1,:)))
end
%% DRUG
numberofExp=2;
LIDDegree=NetParametersExps{numberofExp,1};
DRUGDegree=NetParametersExps{numberofExp,2};
[pMannWhit,h]=ranksum(LIDDegree,DRUGDegree);
[h,pKolmorogov]=kstest2(LIDDegree,DRUGDegree);
NetParam{1}
GEPHIDATA{numberofExp,1}.EXP_ID(1,:)
disp([pKolmorogov;pMannWhit])