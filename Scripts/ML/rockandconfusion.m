% Compute ROC(k) and Confusion Matrix fro several  N model predictions
% 
% Input
%   Ymany:  N model predictions size: observations x models
%   Y:      Ground truth with K classes
% Output
%   ROC:    ROC{1}=[TPR;FPR] vectors fo N models
%   METRICS:Cell kx1 Structure with the following fields:
%       True Positive Rate:     TPR (sensitivity, recall, hit rate)
%       False Positive Rate:    FPR (fall-out)
%       True Negative Rate :    TNR (specificity, selectivity)
%       Positive Predictive Value : PPV (precision)
%       Negative Predictive Value : NPV
%       False Negative Rate :   FNR (miss rate)
%       False Discovery Rate:   FDR
%       False Omission Rate :   FOR
%       F1 score: 2*(PPV*TPR)/(PPV+TPR);
% 
function [ROC,METRICS]=rockandconfusion(Ymany,Y)
[Nobservations,Nmodels]=size(Ymany);
ConditionLabels=unique(Y,'legacy');
Nclasses=numel(ConditionLabels);
Ybin=getbinarylabel(Y); % TARGETS [ground truth]
% MAIN LOOP
ROC=cell(Nclasses,1); % Save FPRs vsTRPs for each Label
METRICS=cell(Nclasses,1); % Save FPRs vsTRPs for each Label
% F1=zeros(ConditionLabels,Nmodels);
for n=1:Nmodels
    Yhat=Ymany(:,n);
    Yhatbin=getbinarylabel(Yhat); % OUTPUTS
    % If classes were unpredicted -> 0 exists
    if ismember(categorical(0),categories(Yhat))
        Yhatbin=Yhatbin(:,unique(Yhat,'legacy')~=categorical(0));
    end
    if ~(size(Ybin,2)==size(Yhatbin,2))
        % Missing predicted class
        MissClasses=setdiff(unique(Y),unique(Yhat));
        for i=1:numel(MissClasses)
            fprintf('>>Unpredicted class: %s\n',char(MissClasses(i)));
            MissIndx(i)=find(ConditionLabels==MissClasses(i));
        end
        OKbin=setdiff(1:Nclasses,MissIndx); % Columns with label
        NewYhatbin=zeros(Nobservations,Nclasses);
        NewYhatbin(:,OKbin)=Yhatbin;
        Yhatbin=NewYhatbin;
    end
    % GET ROC for each LABEL:
    for c=1:Nclasses
        % PORITIVE
        P=sum(Ybin(:,c));
        % NEGATIVE
        N=Nobservations-P;
        %           TRUE        PREDICTED
        % TRUE POSITIVE
        TP=sum((Ybin(:,c)==1).*(Yhatbin(:,c)==1));
        % TRUE NEGATIVE
        TN=sum((Ybin(:,c)==0).*(Yhatbin(:,c)==0));
        % FALSE POSITIVE
        FP=sum((Ybin(:,c)==0).*(Yhatbin(:,c)==1));
        % FALSE NEGATIVE
        FN=sum((Ybin(:,c)==1).*(Yhatbin(:,c)==0));
        % True Positive Rate (sensitivity, recall, hit rate)
        TPR=TP/P;
        % True Negative Rate (specificity, selectivity)
        TNR=TN/N;
        % Positive Predictive Value (precision)
        PPV=TP/(TP+FP);
        % Negative Predictive Value
        NPV=TN/(TN+FN);
        % False Negative Rate (miss rate)
        FNR=FN/P;
        % False Positive Rate (fall-out)
        FPR=FP/N;
        % False Discovery Rate
        FDR=FP/(FP+TP);
        % False Omission Rate
        FOR=FN/(FN+TN);
        % F1 score
        F1=2*(PPV*TPR)/(PPV+TPR);
        % Outputs:
        ROC{c,1}(n,:)=[TPR;FPR];
        METRICS{c,1}.TPR(n)=TPR;
        METRICS{c,1}.FPR(n)=FPR;
        METRICS{c,1}.TNR(n)=TNR;
        METRICS{c,1}.PPV(n)=PPV;
        METRICS{c,1}.NPV(n)=NPV;
        METRICS{c,1}.FNR(n)=FNR;
        METRICS{c,1}.FDR(n)=FDR;
        METRICS{c,1}.FOR(n)=FOR;
        METRICS{c,1}.F1(n)=F1;
        if n==1
            METRICS{c,1}.Label=ConditionLabels(c);
        end
    end
%     [tpr,fpr,thresholds] = roc(Ybin',Yhatbin');
%     [c,cm,ind,per] = confusion(Ybin',Yhatbin');
%       see : >>perfcurve
end