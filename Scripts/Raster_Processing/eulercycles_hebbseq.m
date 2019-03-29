% Function that Look for Euler Cycles at Hebbian Sequence of Neural Ensembles
%
% A cycle is the realization  of all ensembles j,k,..,z from ensemble i-th
% Until ensmeble i-th relizae again
%
% Input:
%   ET: Hebbian Sequence of Neural Ensembles
%   E:  List of ALL Ensembles
% Ouput:
%   TableCycle: Cell :
%       Type of Cycle | Ensemble Sequence: 
% Example:
%       Simple          [1,2,3,1]
%       ...
%       Open            [1,2,3,2,1]
%   TypeCycles: Counter of Type of Cycles: 3 rows vector:
%       TypeCycles(1)-># Simple Cycles
%       TypeCycles(2)-># Closed Cycles
%       TypeCycles(3)-># Open Cycles
% 
% For Instance the Sequence: 1,2,3,2,1,2,3,1,3,1,2 has the sorted cycles:
% Open:     [1,2,3,2,1],[2,3,2,1,2],[3,2,1,2,3]
% Simple:   [1,2,3,1]
% Closed:   [2,3,1,3,1,2]
% NO cycle: [3,1,3,1,2]
% TypeCylces'=[1,1,3]

function [TableCycles,TypeCycles]=eulercycles_hebbseq(ET,E)
% Cycles of Ensembles [REVERBERATION] return to a given ensemble (after activate all ensembles)
TableCycles={};
TypeCycles=zeros(3,1);
tcounter=[]; t=1;
Tremaining=1;
% Combinations of Pairs of Alternances:
aux=1; AllPair=[];
for i=1:numel(E)
    for j=i+1:numel(E)
        AllPair(aux,:)=[E(i),E(j)];
        aux=aux+1;
    end
end
TotalPairs=size(AllPair,1);
auxCyles=1;
while and(~isempty(Tremaining),~isempty(ET))
    ActualEnsemble=ET(t);                   % Get Ensemble by Order
    Cy=find(ET(t+1:end)==ActualEnsemble);   % Following Ensemble Realization Index
    % Check ONly if there were more than 1 Ensemble Realization
    if ~isempty(Cy)     
        auxt=t; % Following Ensembles Index Start
        i=1;    % Index of the Ensemble Realization
        % Eulerian Cycles: all-Ensemble Realizations
        % Search of the Ensembles Between Actual and Folowing
        % Ensembles Realizations and Classify Kind of Cycle:
        while i<=length(Cy)
            % POSSIBLE CYCLE #############################################
            Cycle=ET(auxt:t+Cy(i));
            % Simple: One Realization of each Ensemble :1,2,3,1
            if and(numel(Cycle)==numel(E)+1,numel(unique(Cycle))==numel(E))
                disp(Cycle')
                disp('    \--> Simple Cycle')
                TypeCycles(1)=TypeCycles(1)+1;
                tcounter=[tcounter;t;Cy(1:i)+t];
                auxt=t+Cy(i);
                i=i+1;
                TableCycles{auxCyles,1}='Simple';
                TableCycles{auxCyles,2}=Cycle;
                auxCyles=auxCyles+1;
            elseif numel(unique(Cycle))==numel(E)
                disp(Cycle')
                CycleMat=zeros(max(E));
                for j=1:length(Cycle)-1
                    CycleMat(Cycle(j),Cycle(j+1))=CycleMat(Cycle(j),Cycle(j+1))+1;
                end
                % Pair of Sequences:
                [rowX,colX]=find(CycleMat);
                PairSeq=[rowX,colX];
                CounterPairs=zeros(TotalPairs,1);
                for p=1:numel(rowX)
                    ActualPair=unique(PairSeq(p,:));
                    if numel(ActualPair)>1
                        [~,indexPair]=intersect(AllPair,ActualPair,'rows');
                        CounterPairs(indexPair)=CounterPairs(indexPair)+1;
                    end
                end
                % AllPairs/CouterPairs histogram of pairs of Sequences
                % Check for upper and lower triangles
                % CS=sum(triu(CycleMat)+triu(CycleMat'));
                % If there were sequences from 2 to end
                if sum(CounterPairs>0)==TotalPairs
                    disp('    \--> Closed Cycle')
                    TypeCycles(2)=TypeCycles(2)+1;
                    tcounter=[tcounter;t;Cy(1:i)+t]; % Start & End Index of the Cycle
                    auxt=t+Cy(i); % Next Ensemble Realization Index
                    i=i+1;
                    TableCycles{auxCyles,1}='Closed';
                    TableCycles{auxCyles,2}=Cycle;
                    auxCyles=auxCyles+1;
                else
                    disp('    \--> Open Cycle')
                    TypeCycles(3)=TypeCycles(3)+1;
                    tcounter=[tcounter;t;Cy(1:i)+t];
                    auxt=t+Cy(i);
                    i=i+1;
                    TableCycles{auxCyles,1}='Open';
                    TableCycles{auxCyles,2}=Cycle;
                    auxCyles=auxCyles+1;
                end
            else % Next Ensemble Realization
                disp('>>No Euler Ensembles Activations')
                tcounter=[tcounter;t;Cy(1:i)+t]; % ALL INDEXES CHECKED
                auxt=t+Cy(i); % INcreas Index of the Following Ensemble Realizations
                i=i+1; % Next Ensemble Realization
                % auxt=t;
            end
            tcounter=unique(tcounter);
        end
    else
        Cycle=[];
        tcounter=unique([tcounter;t]);
    end
    Tremaining=setdiff(1:numel(ET)-1,tcounter); % CHECKED INDEXES
    if ~isempty(Tremaining)
        t=Tremaining(1);
    end
end
CyclesLabels={'Simple';'Closed';'Open'};
table(CyclesLabels,TypeCycles)