% Function to extract network parameters from tables
% 
function [TotalTable,DATAnet]=TableNetwork(ExsitFeature,GEPHIDATA,EXPLIST,Names_Conditions,ActualFeature)
TotalTable=table;
DATAnet={};
% ALLpreActive=cell(NE,NC);
fprintf('>>Reading data for table: \n')
n=1;
while sum(ExsitFeature(:))>0
    % No empty Table
    if ExsitFeature(n)
        ActualExp=EXPLIST{n};
        fprintf('> Experiment: %s in ',ActualExp)
        [Rows,Cols]=find(strcmp(ActualExp,EXPLIST));
        preActive=[];
        for k=1:numel(Cols)
            ActualCondition=Names_Conditions{Cols(k)};
            ActualActive=[];
            fprintf('%s,',ActualCondition)
            X=GEPHIDATA{Rows(k),Cols(k)};     % Table
            x=X{:,{ActualFeature{1}}};  % Cell Values
            xdegree= X(:,{'degree'});  % Table Values
            y=X{:,'label'};     % labels bleonging to ensemble
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % Use Nodes with Degree>0 in C_i or peviously Degree>0 in C_i-1
            % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            % Make it Vector: (yeahp, with a loop)
            xnum=zeros(numel(x),1);
            for nx=1:numel(x)
                xnum(nx)=str2double(x{nx});
%                 disp(nx)
                % if y{nx}~='0' && str2double( table2array(xdegree(nx,1)) )>0
                if ~strncmpi(y{nx},'0',numel(y{nx})) && str2double( table2array(xdegree(nx,1)) )>0
                    ActualActive=[ActualActive;nx];
                end
            end
            % Only Previus and Actual Nodes #######################
            OKindex=union(preActive,ActualActive);      % join
            DATAnet{Rows(k),Cols(k)}=xnum(OKindex);
            preActive=OKindex;                  % update
            RowTable=table({ActualCondition},{ActualExp},...
                mean(xnum(OKindex)),mode(xnum(OKindex)),...
                median(xnum(OKindex)),var(xnum(OKindex)),...
                skewness(xnum(OKindex)),kurtosis(xnum(OKindex)));
            % Statistics: mean, variance, mode, median, etc
            TotalTable=[TotalTable;RowTable];
            ExsitFeature(Rows(k),Cols(k))=0;
        end
        fprintf('\n')
    end
    n=n+1;
end
fprintf('>>Data Loaded.\n')