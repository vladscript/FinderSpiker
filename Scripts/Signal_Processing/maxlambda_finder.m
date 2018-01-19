% Function to sparse deconvolve
% It guarantees empty non-zero drivers
% Input
%   XD: Detrended Signal
%   FR: Fluorophore Response 
% Output
%   D:          Driver
%   LAMBDAS:    Sparse Parameter
function [D,LAMBDASS]=maxlambda_finder(XD,FR)

[C,F]=size(XD); % Cells anf Frames
% Initialize Outputs:
D=zeros(C,F);
% XEST=zeros(C,F);
LAMBDASS=zeros(C,1);
% MERRORS=zeros(C,1);
for c=1:C
    x=XD(c,:);
    r=FR(c,:);
    if length(r)>length(x)/2
        disp('Warning: length of response')
        r=r(1:end-length(r)/2);
    end
    d=[];               % Initial Empty Drive
    lambda_pow2 = 10;   % Initial very high 2-power lambda
    % It reduces lambda until it finds a driver signal 
    % Until it starts to fit negative samples
    Positive_Driver=false;
    Negative_Driver=false;
    while ~and(Negative_Driver,Positive_Driver)
        lambda=2^lambda_pow2;
        [d,~,~]=magic_sparse_deconvolution(x,r,lambda);
        d=smooth(d)';%  SMOOTH DRIVER
        if ~isempty( d(d>0) )
            Positive_Driver=true;
        else
            Positive_Driver=false;
        end
        if ~isempty( d(d<0) )
            Negative_Driver=true;
        else
            Negative_Driver=false;
        end
        
        lambda_pow2=lambda_pow2-1;
        disp(' [] [] [] Lambda searcher [] [] []  ')
    end
    % OUTPUtS:
    D(c,:)=d';
    %     XEST(c,:)=x_est';
    %     MERRORS(c)=MSEE;
    % LAMBDASS(c)=lambda^lambda_pow2;
    LAMBDASS(c)=lambda;
end