%% DESCRIPTION *******************************************************
% Script To Make Raster Toy for
% Hierarchichla CLustring of
% Neural Coactivity of
% Multiple & Simultaneous Neural Activity Recording by
% Calcium Imaging
%% SETUP
% Ensmebles
Nensambles=3;       % N ensembles
Nneurons=[3,3,4];   % Neurons at each ensembles
% Shared Neurona Among Ensembles
SharedNeurons=zeros(Nensambles);
SharedNeurons(logical(eye(Nensambles)))=Nneurons;
for n=1:Nensambles-1
    Nrand=round(rand);
    if Nrand
        SharedNeurons(n,n+1)=1;
        SharedNeurons(n+1,n)=1;
    else
        SharedNeurons(n,n+1)=0;
        SharedNeurons(n+1,n)=0;
    end
    % Take Care of Rows are less than Nneurons(n,:)
end
% Hebbian Sequence 
HebbianSeq=[1,2,3,1,2,3,1,2,3,2,1];
%% Build Raster
nens=1:Nensambles;
% Nuerons at each ensemble
% Shared Neurons of Previous Ensembles
TotalNeurons=0;
NeuronsUsed=[];
Ensembles=cell(Nensambles,1);
for n=1:Nensambles
    NSharedafter=sum(SharedNeurons(n,nens(nens>n)),2);
    NSharedbefore=sum(SharedNeurons(n,nens(nens<n)),2);
    Ensembles{n}=TotalNeurons+1-NSharedbefore:TotalNeurons+Nneurons(n)-NSharedbefore;
    NeuronsUsed=[NeuronsUsed,Ensembles{n}];
    TotalNeurons=TotalNeurons+Nneurons(n)-NSharedbefore;
end
NN=max(unique(NeuronsUsed));
Raster=zeros(NN,numel(HebbianSeq));

for n=1:numel(HebbianSeq)
    ActualEnsemble=HebbianSeq(n);
    CoactiveNeurons=Ensembles{ActualEnsemble};
    Raster(CoactiveNeurons,n)=1;
end
%% Permutate
NewIndx=randperm(NN);
Raster=Raster(NewIndx,:);

%% Sample Neurons
%% Analyze
Raster_Analysis=get_bayes_ensembles(Raster);