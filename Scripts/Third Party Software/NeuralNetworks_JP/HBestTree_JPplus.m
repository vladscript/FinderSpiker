% Tree validation
% Measure of how faithfully the tree represents the dissimilarities among
% observations and select the best between different methods: average,
% centroid, complete, median, single, ward and weighted.
%
% HBestTree_JP(Sim, numFig)
%
% Inputs
% Sim = similarity as matrix PxP (P=#peaks).
% numFig = number of the figure to plot.
%
% Outputs
% Validation = validation matrix as 7x2 matrix
%
% ..:: by Jesús E. Pérez-Ortega ::.. Aug-2012
%
% V2.0 output deleted & subplots - modified JP Nov-2012
% Plus Edition in October 2018.
% Plus Edition in January 2019 by github/vladscript.

function recommended = HBestTree_JPplus(Sim)

% method={'average','centroid','complete','median','single','ward','weighted'};
% 'complete' caused MATLAB crashings
% method={'average','ward','weighted'};
method={'average','centroid','median','single','ward','weighted'};

for i=1:numel(method)
    fprintf('> Method Consistency for %s Linkage is:',method{i});
    % Hierarchical Cluster Tree
    Tree = linkage(squareform(1-Sim,'tovector'),method{i});
    
    % Consistency 
    Consistency=inconsistent(Tree);
    Consistency=sum(Consistency(:,4))/...
        (length(Consistency)-length(find(Consistency(:,3)==1)));
    Validation(i)=Consistency;
    fprintf('%3.3f \n',Validation(i));
    
end
[~,Index]=sort(Validation,'descend');
method=method(Index);
recommended=method{1};

% figure(numFig)
% plot(Validation,'-b')
% set(gca,'xtick',1:7,'xticklabel',method)
% title([{['Faithfulness of Hierarchical Tree'];...
%     ['(''' recommended ''' recommended)']}])
% ylabel('Average Consistency')