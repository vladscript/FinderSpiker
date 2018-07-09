% Function to get Raster choosing Method:
% 1: Driver Seignal by Sparse Deconvolution
% 2: Oopsi algorithm
% 3: Derivative method
% Input
% Method: [1,2,3]->{Sparse,oopsi,derivative}
% VARARGIN: according to the method
% Mode 1: Sparse Deconvolution 
%       Driver Signal
%       ActiveNeurons

% Mode 2: Deconvolution
%       Detrended Signal: Only Detected Neurons
%       ActiveNeurons

% Mode 3: Derivative of Cleaned Signal
%       Driver Signal
%       ActiveNeurons
%       Response Function

% Output
% R:        Raster Row vectors of activity
function R=get_raster(method,varargin)
%% Setup Raster Method
switch method
    
    case 1  % Sparse Deconvolutoin (the chido one)
        D=varargin{1};
        ActiveNeurons=varargin{2};
        [TotalCells,F]=size(D);
        R=zeros(TotalCells,F);
        if ~isempty(ActiveNeurons)
            for n=1:numel(ActiveNeurons)
                % [~,Np]=findpeaks(D(n,:)); % way too clean
                R(ActiveNeurons(n),D(ActiveNeurons(n),:)>0)=1;
            end
        end
        
    case 2  % oopsi method
        XD=varargin{1};
        ActiveNeurons=varargin{2}; % vector
        TAUS=varargin{3};
        fs=varargin{4};
        Xest=varargin{5};
        [TotalCells,F]=size(XD);
        R=zeros(TotalCells,F);
        V.dt=1/fs;
        V.smc_iter_max = 1;
        if ~isempty(ActiveNeurons)
            for n=1:numel(ActiveNeurons);
                xd=XD(ActiveNeurons(n),:);
                xe=Xest(ActiveNeurons(n),:);
                P.sig=std(xd-xe);
                tau=TAUS(ActiveNeurons(n),2);
                P.gam   = 1-V.dt/tau;
                y=fast_oopsi(xd,V,P); % CHECK!
                d=zeros(1,F);
                d(y>1)=1;
                R(ActiveNeurons(n),d>0)=1;
            end
        end
        
    case 3  % Derivative Method
        D=varargin{1};
        ActiveNeurons=varargin{2};
        FR=varargin{3};
        [TotalCells,F]=size(D);
        R=zeros(TotalCells,F);
        if ~isempty(ActiveNeurons)
            for n=1:numel(ActiveNeurons)
                d=D(ActiveNeurons(n),:);
                r=FR(ActiveNeurons(n),:);
                x_sparse=sparse_convolution(d,r);
                dx_sparse=diff(x_sparse);
                R(ActiveNeurons(n),dx_sparse>0)=1;
            end
        end
                
    otherwise
        disp('No method'); % not even happening
        
end