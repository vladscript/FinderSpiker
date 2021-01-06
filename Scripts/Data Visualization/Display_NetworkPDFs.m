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
[Ne,Nc]=size(GEPHIDATA);
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
%% Do Magic Plots

for e=1:Ne
    X=GEPHIDATA{e,:};
    figure('Name',sprintf('Experiment: %s ',X.EXP_ID(1,:)),...
        'NumberTitle','off','Position',[20 246 1336 420]);
    MaxK=0;
    MaxParam=0;
    PDF_k=cell(1,Nc);
    PDF_param=cell(1,Nc);
    Nbins_k=zeros(1,Nc);
    Nbins_param=zeros(1,Nc);
    for c=1:Nc
        DataNet=GEPHIDATA{e,c};
        Ndata=numel(DataNet.degree);
        var_k=zeros(Ndata,1);
        netfeat=zeros(Ndata,1);
        % Degree ***************************************************
        for n=1:Ndata
           var_k(n) = str2num( DataNet.degree{n});
           netfeat(n) = str2num( DataNet{:,Ncol}{n} );
        end
        if max(var_k)>0
            k_range=1:max(var_k);
        else
            k_range=unique(var_k);
        end
        MaxK=max([var_k;MaxK]);
        MaxParam=max([netfeat;MaxParam]);
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
end

