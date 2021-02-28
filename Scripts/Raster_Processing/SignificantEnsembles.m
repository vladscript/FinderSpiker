%% Checking Significant Peaks of the Clusteres Raster
% require:
Ensemble_Sorting;

%%
% Input:
cond=2; % Condition
N=4; % Significant level: non random coactivity
% Raster=R_Dyskinesia_Analysis; % Analized Raster
% Raster=R_Clozapine_Analysis; % Analized Raster
Raster=R_Amantadine_Analysis; % Analized Raster

Intervasl=Features_Condition.HebbPositions{cond};
Centers=Features_Condition.HebbCenters{cond};
Hebbs=Features_Condition.HebbianSeq{cond};
% 1) get clusterded frames
Indexes=find(Raster.Peaks.Index); % cluster indexes
CAG=sum(Raster.Data.Data,2);
SingIndex=find(CAG>N);
SingificantEnsembles=[];
SignTimeEns=[];
Start=[];
End=[];
aux=1;
for i=1:numel(Hebbs)
    if ~isempty(intersect( Intervasl(i,1):Intervasl(i,2),SingIndex))
        fprintf('>>Significant Peak @%3.2f [min] Ensemble:%i\n',Centers(i)/fs/60,Hebbs(i));
        SingificantEnsembles(aux,1)=Hebbs(i);
        SignTimeEns(aux,1)=Centers(i)/fs/60;
        Start(aux,1)=Intervasl(i,1);
        End(aux,1)=Intervasl(i,2);
        aux=aux+1;
    end
end
HebbTable=table(SingificantEnsembles,SignTimeEns,Start,End)