% Function to sparse deconvolve by SParSA Algorithm
% It guarantees empty non-zero drivers
% Input
%   XD: Detrended Signal: Row Vectors
%   FR: Fluorophore Response Responses x Length
%   Mode:       arg min{lambda}-> driver:[ +&- or ++ ~=0]
%   no entry    +&- positive and negative driver
%   0           ++  only positive driver
%   1           ~=0 different to zero driver
% Output
%   D:          Driver
%   LAMBDAS:    Sparse Parameter

function [D,LAMBDASS]=maxlambda_finder(XD,FR,varargin)
%% Setup:
% Set Mode to Get Lambda Deconvolution parameter
if isempty(varargin) % [+&- Mode]
    deconmode='posandneg';
elseif varargin{1}==0 % [++ Mode]
    deconmode='onlypos';
elseif varargin{1}==1 % [~=0 Mode]
    deconmode='difftozero';
end
[C,F]=size(XD); % Cells anf Frames
% Initialize Outputs:
D=zeros(C,F);
LAMBDASS=zeros(C,1);
%% Main Loop
for c=1:C
    x=XD(c,:); % Raw signal
    r=FR(c,:); % Response
    if length(r)>length(x)/2
        disp('Warning: length of response')
        r=r(1:end-length(r)/2);
    end
    
    lambda_pow2 = 10;   % Initial very high 2-power lambda
    % It reduces lambda until it finds a driver signal 
    % Until it starts to fit negative samples
    %% First Case: Positive & Negative Driver Signal
    switch deconmode
        case 'posandneg'
            Positive_Driver=false;
            Negative_Driver=false;
            while ~and(Negative_Driver,Positive_Driver)
                if lambda_pow2>0
                    lambda=2^lambda_pow2;
                else
                    lambda=lambda*0.75;
                end
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
                disp(' [*] Lambda searcher Driver Mode: [+] and [-] ')
            end
    %% Second Case: Only Postive Driver Signal
        case 'onlypos'
            d=[]; % Initial Empty Drive
            while isempty(d(d>0))
                if lambda_pow2>0
                    lambda=2^lambda_pow2;
                else
                    lambda=lambda*0.75;
                end
                [d,~,~]=magic_sparse_deconvolution(x,r,lambda);
                % d=smooth(d)';%  SMOOTH DRIVER
                lambda_pow2=lambda_pow2-1;
                disp(' [*] Lambda searcher Driver Mode: [+] ')
            end
    %% Third Case: Just Different to Zero Driver Signal
        case 'difftozero'  
            d=zeros(1,F);               % Initial Zero Drive
            while isempty(d(d~=0))
                if lambda_pow2>0
                    lambda=2^lambda_pow2;
                else
                    lambda=lambda*0.75;
                end
                [d,~,~]=magic_sparse_deconvolution(x,r,lambda);
                % d=smooth(d)';%  SMOOTH DRIVER
                lambda_pow2=lambda_pow2-1;
                disp(' [*] Lambda searcher Driver Mode: [+] OR  [-] ')
            end 
        otherwise
            disp('[No deconvolution done]') % not happening
    end
    
    % OUTPUtS:
    D(c,:)=d';
    %     XEST(c,:)=x_est';
    %     MERRORS(c)=MSEE;
    % LAMBDASS(c)=lambda^lambda_pow2;
    LAMBDASS(c)=lambda;
end