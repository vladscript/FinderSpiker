% Conident predicted labels
% from binary inputs 
function [Yconfident,ClassRatebin]=getrighclassbinary(Y,Yallbin)
Nobser=numel(Y);
labelconditions=unique(Y,'legacy');
Nsim=size(Yallbin{1},2);
% Get right classes:
for k=1:numel(labelconditions)
    Ybin= false(Nobser,1); Yclassbin=Ybin;
    Ybin(Y==labelconditions(k))=true;
    Ymanybin=Yallbin{k};
    % Confident predictions
    ClassRatebin{k}=sum(Ymanybin,2)/Nsim;
    Yclassbin(ClassRatebin{k}==1)=true;
    Yconfident{k}=Yclassbin;
end