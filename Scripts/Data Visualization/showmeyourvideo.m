function showmeyourvideo(X_SPARSE,fs,isSIGNAL,XY,r,W,H,NameOut)
%% Description
% Gets Clean Video Reconstruction from Estimated Fluorescence Signal
% Previously loaded .mat File of Experiment:
% INPUT
%   X_SPARSE=SIGNALSclean{1};   Clean Signals
%   fs:                         Sampling Frequency
%   isSIGNAL = isSIGNALS{i,j};  Index Signals
%   XY:                         Set of Coordinates
%   r:                          Vector of Coordinate Radius
%   W:                          Width   (pixels)  x-axis
%   H:                          Height  (pixels)  y-axis
% OUTPUT
%   Cleansed Video @ NameOut

% EXAMPLE:
% >>load(DatasetTest.mat)
% get_clean_video()
%% Setup *************************************************************
clc; disp('>>Initiazing...')

%% Read IMAGE
% F=size(X_SPARSE,2);
% mov(1:F) = struct('cdata', zeros(H,W,3, 'uint8'),'colormap', []);
% k=floor(F*rand);
% mov(k).cdata = readFrame(VideoObj);
% FRAME=mov(k).cdata;     % Read image


%% ************************************************************************

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

OffFig=20;
figure1=figure('Position', [OffFig, OffFig, W+OffFig, H+OffFig],...
    'Units','pixels');

BLACkGROUND=zeros(H,W,3);
imagesc(BLACkGROUND); colormap gray; hold on;
set(figure1.Children,'Ydir','reverse')
figure1.Children.XLim=[0,W];
figure1.Children.YLim=[0,H];
figure1.Children.Units='pixels';
figure1.Children.Position=[OffFig/2,OffFig/2,W,H];

% Experiment_Condition_Nvideo_1
open(v);
aux_count=1;
for f=1:length(X_SPARSE)
    for i=1:numel(isSIGNAL)
        % XY(i,:);
        % r(i,:);
        xsparse=X_SPARSE(isSIGNAL(i),:);
        xval=xsparse(f);
        if xval<0
            xval=0;
        end
        AmplitutoColor=floor(255*xval/max(clean_x))+1; % GLOBAL
        % AmplitutoColor=floor(255*xval/max(xsparse))+1; % LOCAL
        if isnan(AmplitutoColor); AmplitutoColor=1; end;
        % viscircles(XY(i,:), r(i,:),'Color','b');
        ixy=isSIGNAL(i);
        rectangle('Position',[XY(ixy,1)-r(ixy),XY(ixy,2)-r(ixy),2*r(ixy),2*r(ixy)],...
            'FaceColor',AMP(AmplitutoColor,:),'Curvature',[1 1],'EdgeColor','none');    
        aux_count=aux_count+1;
    end
    %axis([0,W,0,H]); axis off;
    
    % Make Frame after plotting all signals
    % If record=1
    Fmovie=getframe(figure1.Children);
    writeVideo(v,Fmovie);
    % else % preview
%     pause
    % end
    disp(['Recording Video ... ',num2str(100*aux_count/(length(X_SPARSE)*length(isSIGNAL)),5),'%'])
end
close(v); close(figure1);
disp('>>Done.')
