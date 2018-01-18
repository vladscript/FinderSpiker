% Function to detect miss-drivers usually due to artifacts or bleaching
% Input
%   DRIVER:         Matrix of driver signals    [Cells x Frames]
%   X_SPARSE:       Matrix of sparse signals    [Cells x Frames]
% Output
%   Dclean:         Matrix of clean driver      [Cells x Frames]
function [Dclean]=clean_driver(DRIVER,X_SPARSE,XEST)
% Setup
[Cells,~]=size(DRIVER);
Dclean=DRIVER;
for c=1:Cells
    xest=XEST(c,:);         % Estimated signal
    dxest=diff(xest);       % Differential of the estimated
    x_sparse=X_SPARSE(c,:); % Sparse signal
    d=DRIVER(c,:);
    % [~,maxpos]=max(d)
    % Initial missdrive *************************************
    if d(1)>0&&dxest(1)<0
        disp('miss driver @ the beginning');
        stillmissdriver=1;
        Ndet=1;
        while stillmissdriver
            if x_sparse(1+stillmissdriver)>0
                stillmissdriver=1+stillmissdriver;
                % disp('...                   detecting');
                Ndet=Ndet+1;
            else
                stillmissdriver=0;
            end    
        end
        % Find following drivers
        % Find x_sparse turns zero
        Dclean(c,1:Ndet)=0;
    end
    % fINAL missdrive *************************************
    if d(end)>0&&dxest(end)>0
        disp('miss driver @ the beginning');
        stillmissdriver=1;
        Ndet=1;
        while stillmissdriver
            if x_sparse(1+stillmissdriver)>0
                stillmissdriver=1+stillmissdriver;
                % disp('...                   detecting');
                Ndet=Ndet+1;
            else
                stillmissdriver=0;
            end    
        end
        % Find following drivers
        % Find x_sparse turns zero
        Dclean(c,1:Ndet)=0;
    end
end