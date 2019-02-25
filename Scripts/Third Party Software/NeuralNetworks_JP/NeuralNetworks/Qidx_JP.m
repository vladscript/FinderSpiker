% Modularity
% Get modularity.
%
% [ClustIdx g] = Qidx_JP(Tree, Data, Adjacent, numFig)
%
% Inputs
% Tree = hierarchical cluster tree
% Adjacent = correlation as a square matrix
% numFig = number of the figure to plot
%
% Output
% ClustIdx = clustering indexes from 2 to 10 groups
% g = best number of groups according to the index selected
%
% ..:: by Jesús E. Pérez-Ortega ::.. April-2013

function [ClustIdx g] = Qidx_JP(Tree, Data, Adjacent, numFig)

Adjacent=double(Adjacent>0);
total_edges=sum(sum(Adjacent))/2;
N=length(Adjacent);

ClustIdx=zeros(10,1);
for i=2:10
    T = cluster(Tree,'maxclust',i);
    g=max(T);
    
    %{
    G=zeros(1,N);
    for j=1:g
        idx=find(T==k);
        g1=find((sum(Data(idx,:),1)/length(idx))>.7);
        G(g1,j)=j;
    end
    %}
    
    
    % Edges between communities
    E=zeros(i);
    for j=1:g
        idx=find(T==j);
        g1=find(sum(Data(idx,:),1)/length(idx));
        n_g1(j)=length(g1);
        for k=1:g
            idx=find(T==k);
            g2=find(sum(Data(idx,:),1)/length(idx));
            %
            if j~=k
                g2=setdiff(g1,g2);
            end
            %}
            edges=sum(sum(Adjacent(g1,g2)))/2;
            E(j,k)=edges;
            E(k,j)=edges;
        end
    end
    
    Q=0;
    for j=1:g
        Dint=2*E(j,j)/(n_g1(j)*(n_g1(j)-1));
        ext=setdiff(1:g,j);
        Dext=sum(E(j,ext))/(n_g1(j)*(N-n_g1(j)));
        Qj(j)=Dint-Dext;
    end
    Q=mean(Qj);%/std(Qj);
    %{
    E=E/total_edges;
    Q=trace(E)-sum(sum(E^2));
    %}
    ClustIdx(i)=Q;
end


figure(numFig)
plot(ClustIdx,'-b');
hold on
[gIdx g]=max(ClustIdx);
plot(g,gIdx,'*r')
title(['Mean Modularity (' num2str(g) ' groups recommended)'])