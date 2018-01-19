%% Function to load video of fluorescence .avi
% Input
%   FileName: Name of the video
%   PathName: Directory of the video
% Output
%   mov: Video, Structure with pixel's fluorescence (grayscale)
function [mov]=Video_Load(FileName,PathName)
%% Read Movie
disp('****Loading video***')
% [FileName,PathName] = uigetfile('*.avi','Select the videos.');
VideoObj = VideoReader([PathName FileName]);
F = VideoObj.Duration*VideoObj.FrameRate; % Frames
W = VideoObj.Width;                       % Weight
H = VideoObj.Height;                      % Height
% FWH=[F,W,H];
mov(1:F) = struct('cdata', zeros(H,W,3, 'uint8'),'colormap', []);
k = 1;
% h=waitbar(0,['Loading Video',FileName]);
while hasFrame(VideoObj)
    disp(['Video: ',FileName(1:end-4),' frame: ',num2str(k),'/',num2str(F)])
    mov(k).cdata = readFrame(VideoObj);
    k = k+1;
%     waitbar(k/F)
end
% close(h)
