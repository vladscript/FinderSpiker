% Funtiio to accumulate Rate of Activity of Each Cell
% For each Experiment
% Input
[~,NC]=size(R_Condition);
% RoA_ALL=cell(1,NC);
% RoA_POS=cell(1,NC);
% RoA_NEG=cell(1,NC);
for c=1: NC
    disp(c)
    R_ALL=R_Condition{c};
%     R_POS=R_merged{c};
%     R_NEG=R_nomerged{c};
    [N_ALL,F_ALL]=size(R_ALL);
%     [N_POS,F_POS]=size(R_POS);
%     [N_NEG,F_NEG]=size(R_NEG);
%     AllActive=[MetaDataColocaliation.PositiveCells;MetaDataColocaliation.NegativeCells];
%     PosActive=MetaDataColocaliation.PositiveCells;
%     NegActive=MetaDataColocaliation.NegativeCells;
    % Cummulate Stuff:
    AllActive=find(sum(R_ALL,2));
    RoA_ALL{c}=[RoA_ALL{c}();sum(R_ALL(AllActive,:),2)./F_ALL];
%     RoA_POS{c}=[RoA_POS{c};sum(R_merged{c},2)./F_POS];
%     RoA_NEG{c}=[RoA_NEG{c};sum(R_nomerged{c},2)./F_NEG];
end
disp('end')