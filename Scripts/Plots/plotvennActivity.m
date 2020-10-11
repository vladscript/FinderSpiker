%% Plot Venn Diagram
%  Input
%   ActivityGroups
% Output
%   Plot
function plotvennActivity(ActivityGroups)
% 3 main SEPARATED groups:
KeepDoPlot=true;
if isempty(intersect(intersect( ActivityGroups.Nfac, ActivityGroups.Nunc),...
ActivityGroups.Ndep))
    fprintf('>>OK')
    AllCells=union(union( ActivityGroups.Nfac, ActivityGroups.Nunc),ActivityGroups.Ndep);
else
    fprintf('Something Wrong.\n');
    KeepDoPlot=false;
end

% Subgroups of

if sum(ismember(ActivityGroups.Nact,ActivityGroups.Nfac))==numel(ActivityGroups.Nact);
    fprintf('>>OK')
else
    fprintf('Something Wrong.\n');
    KeepDoPlot=false;
end

% ActivityGroups.Ninh % is IN:
% ActivityGroups.Ndep

if sum(ismember(ActivityGroups.Ninh,ActivityGroups.Ndep))==numel(ActivityGroups.Ninh);
    fprintf('>>OK')
else
    fprintf('Something Wrong.\n');
    KeepDoPlot=fasle;
end

if KeepDoPlot
    %% Areas & Intersections
    TotalArea=300;
    Afac=numel(ActivityGroups.Nfac)/numel(AllCells)*TotalArea;
    Aunc=numel(ActivityGroups.Nunc)/numel(AllCells)*TotalArea;
    Adep=numel(ActivityGroups.Ndep)/numel(AllCells)*TotalArea;
    Inter=[0,0,0,0];
    figure %, axis equal, axis off
    venn([Afac,Aunc,Adep],Inter,'FaceColor',{'g','b','r'},'FaceAlpha',{1,1,1},'EdgeColor','black')
    legend('Facilitated','Unchanged','Depressed');
    axis image
    % Subgroups of Fac & Dep
    Aact=numel(ActivityGroups.Nact)/numel(ActivityGroups.Nfac)*Afac;
    Ainh=numel(ActivityGroups.Ninh)/numel(ActivityGroups.Ndep)*Afac;
    figure %, axis equal, axis off
    venn([Afac,Aact],Aact,'FaceColor',{[0,1,0],[0.5,1,0.5]},'FaceAlpha',{1,0.6},'EdgeColor','black')
    legend('Facilitated','Activated');
    axis image
    figure% , axis equal, axis off
    venn([Adep,Ainh],Ainh,'FaceColor',{[1,0,0],[1,0.5,0.5]},'FaceAlpha',{1,0.6},'EdgeColor','black')
    legend('Depressed','Inhibited');
    axis image
end