% Function that Look for Euler Cycles at
% Hebbian Sequence of Neural Ensembles
% Input_
%   ET: Hebbian Sequence of Neural Ensembles
% Ouput_
%   TypeCycles: Counter of Type of Cycles:
% TypeCycles(1)->Simple
% TypeCycles(2)->Closed
% TypeCycles(3)->Open
function TypeCycles=eulercycles_hebbseq(ET,E)
% Cycles of Ensembles [REVERBERATION] return to a given ensemble (after activate all ensembles)
TypeCycles=zeros(3,1);
tcounter=[]; t=1;
Tremaining=1;
while and(~isempty(Tremaining),~isempty(ET))
    ActualEnsemble=ET(t);                   % Get Ensemeble
    Cy=find(ET(t+1:end)==ActualEnsemble);   % Ensemble Realization Frames
    % Check ONly if there were more than 1 Ensemble Realization
    if ~isempty(Cy)     
        auxt=t;
        i=1;
        % Eulerian Cycles: all-Ensemble Realizations
        % Search of the Ensembles Between Actual and Folowing
        % Ensembles Realizations and Classify Kind of Cycle:
        while i<=length(Cy)
            Cycle=ET(auxt:t+Cy(i));
            % Simple: One Realization of each Ensemble :1,2,3,1
            if and(numel(Cycle)==numel(E)+1,numel(unique(Cycle))==numel(E))
                disp(Cycle')
                disp('Simple')
                TypeCycles(1)=TypeCycles(1)+1;
                tcounter=[tcounter;t;Cy(1:i)+t];
                auxt=t+Cy(i);
                i=i+1;
            elseif numel(unique(Cycle))==numel(E)
                disp(Cycle')
                CycleMat=zeros(max(E));
                for j=1:length(Cycle)-1
                    CycleMat(Cycle(j),Cycle(j+1))=CycleMat(Cycle(j),Cycle(j+1))+1;
                end
                % Check for upper and lower triangles
                CS=sum(triu(CycleMat)+triu(CycleMat'));
                % If there were sequences from 2 to end
                if sum(CS(2:end)>0)==numel(E)-1
                    disp('Closed')
                    TypeCycles(2)=TypeCycles(2)+1;
                    tcounter=[tcounter;t;Cy(1:i)+t];
                    auxt=t+Cy(i);
                    i=i+1;
                else
                    disp('Open')
                    TypeCycles(3)=TypeCycles(3)+1;
                    tcounter=[tcounter;t;Cy(1:i)+t];
                    auxt=t+Cy(i);
                    i=i+1;
                end
            else % Next Ensemble Realization
                disp('No Euler Activations')
                tcounter=[tcounter;t;Cy(1:i)+t];
                i=i+1;
                % auxt=t;
            end
            tcounter=unique(tcounter);
        end
    else
        Cycle=[];
        tcounter=unique([tcounter;t]);
    end
    Tremaining=setdiff(1:numel(ET)-1,tcounter);
    if ~isempty(Tremaining)
        t=Tremaining(1);
    end
end