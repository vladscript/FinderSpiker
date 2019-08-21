%% DESCRIPTION *******************************************************
% Script To Make Raster Toy for
% Hierarchichal Clustring of
% Neural Coactivity for
% Multiple & Simultaneous Neural Activity by
% Calcium Imaging Recording

%% NEST TASKS
% Add Factor of CLuster: 0.5,1-> When Creatiing Ensembles

%% SETUP
clc; clear;
disp('>>Building Raster...')
% Ensmebles
Nensambles=5;       % N ensembles
Nneurons=[3,6,7,2,5];   % Neurons at each ensembles

% Hebbian Sequence 
HebbianSeq=[1,2,3,4,5,1,2,3,4,5,5,4,3,2,1];

% Noisy frames
Frames=numel(HebbianSeq);
LevelRandFrame=0.05; % [0,1]

% Random Neurons
LevelRandNeu=0.05;   % [0,1]

% Ratio of random activity
LevelNoise=0.1; % % Active Neurons or Frames at random Activity


% Shared Neurons Among Ensembles (Random)
SharedNeurons=zeros(Nensambles);
SharedNeurons(logical(eye(Nensambles)))=Nneurons;
for n=1:Nensambles-1
    Nrand=round(rand);
    MaxSharinNeu=floor(Nneurons(n)/2);
    Neuron2Share=0:MaxSharinNeu;
    Neuron2ShareIndx=randperm(MaxSharinNeu+1,1);
    if Nrand
        SharedNeurons(n,n+1)=Neuron2ShareIndx;
        SharedNeurons(n+1,n)=Neuron2ShareIndx;
    else
        SharedNeurons(n,n+1)=0;
        SharedNeurons(n+1,n)=0;
    end
    % Take Care of Rows are less than Nneurons(n,:)
end


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
% Set dots at Raster
for n=1:numel(HebbianSeq)
    ActualEnsemble=HebbianSeq(n);
    CoactiveNeurons=Ensembles{ActualEnsemble};
    Raster(CoactiveNeurons,n)=1;
end

%% Add Random Activity
RasterRandom=Raster;
NrandFrames=floor(LevelRandFrame*Frames);
NrandNeuron=round(LevelRandNeu*NN);
% Random Columns
Rrandcol=zeros(NN,NrandFrames);
Nrandneurons=floor(LevelNoise*NN);
for n=1:NrandFrames
    IndexRandNeuron=randperm(NN,Nrandneurons);
    Rrandcol(IndexRandNeuron,n)=1;
    % Add Random Column to the Raster
    AfterFrame=randperm(Frames,1);
    RasterRandom=[RasterRandom(:,1:AfterFrame-1),Rrandcol(:,n),RasterRandom(:,AfterFrame:end)];
end
% Random Rows
Rrandrow=zeros(NrandNeuron,Frames+NrandFrames);
NrandFramed=round(LevelNoise*(Frames+NrandFrames));
for n=1:NrandNeuron
    IndexRandFrames=randperm(Frames+NrandFrames,NrandFramed);
    Rrandrow(n,IndexRandFrames)=1;
end
RasterRandom=[RasterRandom;Rrandrow];




%% Plot DATA
disp('>>Built Raster.')
disp('Resume Ensembles: ')
tabulate(HebbianSeq);
disp('Press ENTER to Show Rasters...')
pause;
%% PLOT & ANALYZE CLEAN DATA
Plot_Raster_Ensembles(Raster);
disp('Example of Active Neurons')
NewIndx=randperm(NN);
Raster=Raster(NewIndx,:);
disp('Permutted Active Neurons')
Plot_Raster_Ensembles(Raster,1,1,NewIndx);
disp('>> Press any key To Analyze');
pause;
Raster_Analysis=get_bayes_ensembles(Raster);
disp('>>Show Data + Random Activity');
pause;

%% PLOT & ANALYZE CLEAN DATA + RANDOM ACTIVITY

Plot_Raster_Ensembles(RasterRandom);
disp('Example of Active Neurons')
NewIndx=randperm(NN+NrandNeuron);
RasterRandom=RasterRandom(NewIndx,:);
disp('Permutted Active Neurons')
Plot_Raster_Ensembles(RasterRandom,1,1,NewIndx);
disp('>> Press any key To Analyze');
pause;
RasterRandom_Analysis=get_bayes_ensembles(RasterRandom);
%% END
fprintf('\n               In the Algorithm We Trust\n\n')