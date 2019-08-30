% Function that gets the indexes to read variables from workspace
% with a specific sufix -Filterby-
% Input
%   WorkSpaceVariables: list of variable names (who function)
%   Filterby: Sufix of the Variable
%   NC: N variables to read
% Output
%   index_var   index of the workspace variable
%   index_CHECK check if it was picked or ignored (0)
function [index_var,WorkspaceVariables]=get_var_index(WorkspaceVariables,Filterby,NC)
% index_var=zeros(1,NC);
index_var=cell(NC,1);
index_CHECK=index_var;
[Nvar,~]=size(WorkspaceVariables);
Variable_Index_Filtered=[];
for nv=1:Nvar
    NamVar=WorkspaceVariables{nv};
    SameWords=ismember(Filterby,NamVar);
    if sum(SameWords)==length(Filterby)
        Variable_Index_Filtered=[Variable_Index_Filtered,nv];
    end
end
WorkspaceVariables=WorkspaceVariables(Variable_Index_Filtered);
% Pick Variables Analysis to Plot
for c=1:NC
    [index_var{c},index_CHECK{c}] = listdlg('PromptString','Select a _Analysis Variable:',...
                    'SelectionMode','single',...
                    'ListString',WorkspaceVariables);
end