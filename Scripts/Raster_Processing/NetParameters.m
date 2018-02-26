% Function to Calculate Network Parameters
% Input
%   Raster:         Activity Matrix: Cells x Frames
%   th_frames:      Threshold of Links - Frames
% Ouput
%   A:              Adjacency Matrix: Cells x Cells
%   kdegree:        K Links for each Cell
%   pk:             links pdf
%   ck:             clustering coefficient
%   NetFeatures:    PowerLaw Parameters & Fitting: param1*k~(-param2)| r^2: param3
function [A,kdegree,pk,binK,ck,NetFeatures]=NetParameters(Rx,th_frames)
    % Adjacency Matrix Data:
    A=GetAdjacencyMatrix(Rx);           % A(i,j):N frames when Cells i & j Fired Together
    Abin=zeros(size(A));                % Initialize Binary  Adjacency Matrix                
    Abin(A>=th_frames)=1;               % Binary  Adjacency Matrix
    kdegree=sum(Abin);                  % k-degree: number of links for each ccell
    % [pk,kbin]=ksdensity(kdegree,'support','positive','function','survivor');

    % Network Parameters
    % k Histogram as proxy to p(k)~a*k^-gamma
    [CountsK,bink]=histcounts(kdegree);
    pk=CountsK/sum(CountsK);
    binK=(bink(1:end-1)+bink(2:end))/2;
    % Peak of pdf:
    [~,maxpk]=max(pk);
    if maxpk==min(binK)
        [fitresult, gof] = powerlawfit(binK(CountsK>0),pk(CountsK>0));
    else
        [fitresult, gof] = powerlawfit(binK(binK>maxpk-1),pk(binK>maxpk-1));
    end
    gamma_k=fitresult.b;
    A_pk=fitresult.a;
    rsquare_pk=gof.rsquare;
    % Clustering Coefficient: c(k)~b*k^-alpha
    C=clustering_coef_bu(Abin); % Degree of Neighbours Interoconections
    kbin=unique(kdegree);
    ck=zeros(size(kbin));
    for i=1:length(kbin)
        currentneruons= kdegree==kbin(i);
        ck(i)=mean(C(currentneruons));
    end
    [fitresult, gof] = powerlawfit(kbin(kbin>1),ck(kbin>1));
    alpha_k=fitresult.b;
    B_ck=fitresult.a;
    rsquare_ck=gof.rsquare;
    % Parameters Output: Expoenent parameter>0
    NetFeatures=[A_pk,-gamma_k,rsquare_pk;B_ck,-alpha_k,rsquare_ck];
end