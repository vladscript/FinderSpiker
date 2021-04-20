% Generate Table of Features for Merged and NO Merged Cells
% Input
%   Raster w/dimension: Cells x Frames
%   Type of Cells: Merged or No Merged
% Output
%   CSV: File with all the Features
function DATApdf=gen_feat_table_merged(R,fs,typecells,Experiment,Names_Conditions)
%% GET FEATURES
% [~,~]=size(R);
[Descriptive,AI_1,AI_2]=get_general_features(R);
[ISIbin,ISIp,Lengthp,Lengthbin,StatsFeatures]=get_iti_pdf(R,fs);

RAN=Descriptive.RAN;
DurAc=Descriptive.DurAc;    
CAG=Descriptive.CAG;
RoA=Descriptive.RoA;

DATApdf.isibin=ISIbin;
DATApdf.isip=ISIp;
DATApdf.LTp=Lengthp;
DATApdf.LTbin=Lengthbin;
DATApdf.cag=CAG;
DATApdf.roa=RoA;
% Stats...
RoAstats=[mean(RoA),mode(RoA),var(RoA),skewness(RoA),kurtosis(RoA)];
CAGstats=[mean(CAG),mode(CAG),var(CAG),skewness(CAG),kurtosis(CAG)];
%% CREATE TABLE
HeadersFeatures={'RateNeurons','ActivityTimeFraction','MeanActivity','EffectiveActivity',...
        'ISImean','ISImode','ISIvar','ISIskew','ISIkurt',...
        'Lengthmean','Lengthmode','Lengthvar','Lengthskew','Lengthkurt',...
        'CAGmean','CAGmode','CAGvar','CAGskew','CAGkurt',...
        'RoAmean','RoAmode','RoAvar','RoAskew','RoAkurt'};
    
Trasterfeatures=table(RAN,DurAc,AI_1,AI_2,...
        StatsFeatures(1),StatsFeatures(2),StatsFeatures(3),StatsFeatures(4),StatsFeatures(5),...
        StatsFeatures(6),StatsFeatures(7),StatsFeatures(8),StatsFeatures(9),StatsFeatures(10),...
        CAGstats(1),CAGstats(2),CAGstats(3),CAGstats(4),CAGstats(5),...
        RoAstats(1),RoAstats(2),RoAstats(3),RoAstats(4),RoAstats(5),...
        'VariableNames',HeadersFeatures);
    
disp(Trasterfeatures);

%% SAVE CSV FILE

FileDirSave=pwd;
slashes=find(FileDirSave=='\');
FileDirSave=FileDirSave(1:slashes(end));

NameDir='Raster Features\';
if isdir([FileDirSave,'\Raster Features'])
    writetable(Trasterfeatures,[FileDirSave,NameDir,Experiment(2:end),'_',Names_Conditions,typecells,'_Features.csv'],...
        'Delimiter',',','QuoteStrings',true);
    disp(['Saved Merging Features: ',Experiment,'-',Names_Conditions])
else % Create Directory
    disp('Directory >Raster Features< created')
    mkdir([FileDirSave,NameDir]);
    writetable(Trasterfeatures,[FileDirSave,NameDir,Experiment(2:end),'_',Names_Conditions,typecells,'_Features.csv'],...
        'Delimiter',',','QuoteStrings',true);
    disp('Resume Tables Directory Created');
    disp(['Saved Merging Features: ',Experiment,'-',Names_Conditions])
end
%% END OF THE WORLD