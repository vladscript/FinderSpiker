% Count if labels were missclassified in N models
% Input
%   Ymany:      Predicted Labels Matrix, size=Nobservation x Nmodels 
%   Y:          Ground Truth Labels
% Output
%   Yclass:     Vector with well-classified in ALL models 
%   ClassRate:  Rate of well classified
function [Yclass,ClassRate]=getrightclass(Ymany,Y)
% Setup
UniqueLabels=unique(Y,'legacy');
Nlabels=numel(UniqueLabels);
[Nobs,Nmods]=size(Ymany);
ClassRate=zeros(Nobs,Nlabels);
Yclass=categorical(zeros(Nobs,1));
% Main loop
for n=1:Nobs
    Nprediction=(Ymany(n,:));
    for c=1:Nlabels
        ClassRate(n,c)=sum(Nprediction==UniqueLabels(c))/Nmods;
    end
end
% Find confusing predictions:
Confused=find(max(ClassRate')<1);
Secure=find(max(ClassRate')==1);
Yclass(Secure)=Ymany(Secure,end);

