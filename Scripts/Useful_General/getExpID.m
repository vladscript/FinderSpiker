function Experiment=getExpID(PathName)
%% Get the Experiment ID:

slashes=find(PathName=='\');
if ~isempty(slashes)
    Experiment=PathName(slashes(end-1)+1:slashes(end)-1); % Experiment ID
else
    Experiment=PathName(1:find(PathName=='.',1)-1);
end
Experiment(Experiment==' ')='_'; % REPLACE SpaceS with '_'
ExpIDDef{1}=Experiment;
% Confirm ID of the Experiment
ExpInput= inputdlg('Experiment ID: ','Confirm unique ID:/Press Canel to use default', [1 75],ExpIDDef);
if ~isempty(ExpInput)
    Experiment=ExpInput{1};
    Experiment(Experiment==' ')='_'; % REPLACE SpaceS with '_'
end
fprintf('>>Experiment ID: %s\n',Experiment)    
