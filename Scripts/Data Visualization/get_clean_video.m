function get_clean_video(X_SPARSE,fs,isSIGNAL,XY,r,W,H,NameOut)
%% Description
% Gets Clean Video Reconstruction from Estimated Fluorescence Signal
% Previously loaded .mat File of Experiment:
% INPUT
%   i:condition; j:#video from Experiment: eg.:RASTER{j,i};
%   X_SPARSE=SIGNALSclean{1};   Clean Signals
%   isSIGNAL = isSIGNALS{i,j};  Index Signals
%   XY:                         Set of Coordinates
%   r:                          Vecor of Coordinate Radius
%   W:                          Width   (pixels)  x-axis
%   H:                          Height  (pixels)  y-axis
% OUTPUT
%   Cleansed Video @ ....

% EXAMPLE:
% >>load(DatasetTest.mat)
% get_clean_video()
%% Setup *************************************************************
clc; disp('>>Initiazing...')

%% Read IMAGE
F=size(X_SPARSE,2);
mov(1:F) = struct('cdata', zeros(H,W,3, 'uint8'),'colormap', []);
% k=floor(F*rand);
% mov(k).cdata = readFrame(VideoObj);
% FRAME=mov(k).cdata;     % Read image


%% ************************************************************************
% *************************************************************************
% Cleaned Video ***********************************************************
% *************************************************************************
% *************************************************************************

% Previosly read X_SPARSE
% isSIGNAL=isSIGNALS{1};
x_sparse_signal=X_SPARSE(isSIGNAL,:);
clean_x=x_sparse_signal(:);
% Three Rule Magic: ********************
% x<=0          ->     0
% max(clean_x)  ->   2^Bits=2^8=256
% xi            ->   256*xi/max(clean_x)

v = VideoWriter([NameOut,'.avi']);
v.FrameRate=fs;
% v.Height=H;
% v.Width=W;

AMP=gray(256);                  % Colormap

figure1=figure('Position', [0, 0, W+146, H+116]);
BLACkGROUND=zeros(H,W,3);
imagesc(BLACkGROUND); colormap gray; hold on;
set(gca,'Ydir','reverse')
% Experiment_Condition_Nvideo_1
open(v);
aux_count=1;
for f=1:length(X_SPARSE)
    for i=1:length(isSIGNAL)
        % XY(i,:);
        % r(i,:);
        xsparse=X_SPARSE(isSIGNAL(i),:);
        xval=xsparse(f);
        if xval<0
            xval=0;
        end
        % AmplitutoColor=floor(255*xval/max(clean_x))+1; % GLOBAL
        AmplitutoColor=floor(255*xval/max(xsparse))+1; % LOCAL
        % viscircles(XY(i,:), r(i,:),'Color','b');
        rectangle('Position',[XY(i,1)-r(i),XY(i,2)-r(i),2*r(i),2*r(i)],...
            'FaceColor',AMP(AmplitutoColor,:),'Curvature',[1 1],'EdgeColor','none');    
        aux_count=aux_count+1;
    end
    axis([0,W,0,H]); axis off;
    
    % Make Frame after plotting all signals
    % If record=1
    Fmovie=getframe(figure1);
    writeVideo(v,Fmovie);
    % else % preview
%     pause
    % end
    disp(['Recording Video ... ',num2str(100*aux_count/(length(X_SPARSE)*length(isSIGNAL)),5),'%'])
end
close(v)
disp('>>Done.')
