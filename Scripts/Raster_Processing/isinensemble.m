% function to idnetify if a neural ensembles is a subset 
% from another neural ensemble
% Input:
%   NeuroVectors:   Matrix where 1 is for active Neuron Index
% Ouput
%   EnsemblesIndx:  Columns of Ensembles which belong to Other: 
%                   [i,j]: ith ensemble is member of the jth ensemble
function EnsemblesIndx=isinensemble(NeuroVectors)
% Initialize
% isSubset=false;
EnsemblesIndx=[];
[~,Nensambles]=size(NeuroVectors);
fprintf('Analysing %i Ensembles ',Nensambles);
% Main Loop
for n=1:Nensambles
    ActualEnsemble=NeuroVectors(:,n);
    ElseEnsemblesIndex=setdiff([1:Nensambles],n);
    for m=1:numel(ElseEnsemblesIndex)
        TestEnsembles=NeuroVectors(:,ElseEnsemblesIndex(m));
        % Ask if actual ensembles belongs to test (else) ensemble(s) 
        % if sum(ismember(find(ActualEnsemble),find(TestEnsembles)))==numel(find(TestEnsembles))
        if sum(ismember(find(TestEnsembles),find(ActualEnsemble)))==numel(find(ActualEnsemble))
            %isSubset=true;
            EnsemblesIndx=[EnsemblesIndx;n,ElseEnsemblesIndex(m)];
        end
        fprintf('.');
    end
end
fprintf('\n');