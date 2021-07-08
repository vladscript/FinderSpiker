% Function to detect empty tables and empty columns of netowrok fetures
% 
% 
function [ExsitFeature,EXPLIST]=emptycolnettables(GEPHIDATA,Names_Conditions,ActualFeature)
% Init Output
[NE,NC]=size(GEPHIDATA); % Number of Experiments & Conditions
ExsitFeature=zeros(NE,NC);
EXPLIST=cell(NE,NC);
for e=1:NE
    for c=1:NC
        X=GEPHIDATA{e,c};
        if ~isempty(X)
            EXPID=X{1,2};
            CurrentFeatures=X.Properties.VariableNames(6:end);
            ExsitFeature(e,c)=1;
        else
            fprintf('>>Missing Experiment @ %s\n',Names_Conditions{c})
            EXPID='no_exp';
            CurrentFeatures={'empty'};
        end
        EXPLIST{e,c}=EXPID;

        ColumndIndx=strmatch(ActualFeature,CurrentFeatures);
        if ~isempty(ColumndIndx)
            ExsitFeature(e,c)=1;
            fprintf('*')
        else
            fprintf('>>Missing Network parameter: %s @ EXP: %s\n',ActualFeature{1},EXPLIST{e,c});
        end
    end
    fprintf('\n')
end