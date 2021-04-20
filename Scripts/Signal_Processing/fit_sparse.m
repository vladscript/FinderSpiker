% Function that fits deconvolution
% Input
%   X_SPARSE
%   Xest
%   XDupdate
%   FR
%   LAMBDASS
%   DRIVER
%   ActiveNeurons
% Output
%   LAMBDASS
%   X_SPARSE
%   DRIVER
function [LAMBDASS,X_SPARSE,DRIVER]=fit_sparse(X_SPARSE,Xest,XDupdate,LAMBDASS,FR,DRIVER,ActiveNeurons)
% AllStdNoises=std(XDupdate'-Xest');
for c=1:numel(ActiveNeurons)
    signdx=ActiveNeurons(c);
    x=XDupdate(signdx,:);
    xdenoised=Xest(signdx,:);
    xsparse=X_SPARSE(signdx,:);
    r=FR(signdx,:);
    lamb=LAMBDASS(signdx,1);
    pred=DRIVER(signdx,:);
    % plot(xdenoised); hold on;
    % plot(xsparse); hold off;
    % big lambdas
    if lamb>1 && std(x-xsparse)>std(x-xdenoised) 
        fprintf('>>Set Unitary Lambda:\n');
        oklambda=1;
        [d,~,~]=magic_sparse_deconvolution(x,r,oklambda);
        xsk=sparse_convolution(d,r);
        % Update Data
        X_SPARSE(signdx,:)=xsk';
        LAMBDASS(signdx,1)=oklambda;
        % noisden=x-xdenoised;
        % noisparse=x-xsk';
        % std(noisparse)
        fprintf('\n>>Set Unitary Lambda: done For Signal %i.\n',signdx);
    elseif std(x-xsparse)>std(x-xdenoised) 
        fprintf('>>Diminish Lambda:\n');
        oklambda=lamb/2;
        [d,~,~]=magic_sparse_deconvolution(x,r,oklambda);
        xsk=sparse_convolution(d,r);
        % Update Data
        X_SPARSE(signdx,:)=xsk';
        LAMBDASS(signdx,1)=oklambda;
        fprintf('\n>>Diminish Lambda: done For Signal %i.\n',signdx);
    else
        fprintf('\n>>Nice Lambda For Signal %i.\n',signdx);
        d=pred;
    end
    % CHECK AND CLEAN DRIVER
    fprintf('\n>>Checking Driver Signal %i.\n',signdx);
    [difx,~,~,~,~,~]=analyze_driver_signal(d,r,x,xdenoised);
    fprintf('\n>>Checked Driver Signal %i.\n',signdx);
    DRIVER(signdx,:)=difx;
end
fprintf('\n****************************************\n\n')