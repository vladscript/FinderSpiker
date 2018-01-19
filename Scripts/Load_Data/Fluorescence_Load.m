%% Function to load fluorescence signals from videos
%  using the ring acquisition method
% Input
%   mov: video of grayscale pixels
%   XY: Coordinates of active signals
%   r: vector of radius of the coordinates
% Output
%   FluoSignals: Matrix of Fluorescence Signals
function [FS]=Fluorescence_Load(mov,XY,r)
%% Start
FS=[];                          % Fluorescence Signals
F=length(mov);                  % Number of frames
[H,W]=size(mov(1).cdata);       % Size of the frames
[NS,~]=size(XY);                % Number of Signals
%     figure(1);
%     imagesc(mov(floor(F*rand)).cdata);
%     hold on; scatter(XY(:,1),XY(:,2),pi*r.^2,'LineWidth',1)
%     title('Loading Flurorescence'); axis off
%     disp('****Loading Signals***')
for n=1:NS      % for every coordinate
    % Circle's Area Pixels ***************************************
    Mx=XY(n,1)-(r):XY(n,1)+(r); % range in x of square
    My=XY(n,2)-(r):XY(n,2)+(r); % range in y of square
    Meshxy_Circle=[]; aux1=1;
    for i=1:length(Mx)
        % chech if it's in image's limits Xaxis
        if Mx(i)>0 && Mx(i)<=W 
            for j=1: length(My)
                % chech if it's in image's limits Yaxis
                if My(j)>0 && My(j)<=H 
                    % check if it's in circle
                    if (Mx(i)-XY(n,1))^2+(My(j)-XY(n,2))^2<=r(n)^2
                        Meshxy_Circle(aux1,:)=[Mx(i),My(j)];
                        aux1=aux1+1;
                    end
                end
            end
        end
    end
    % Fluorescence ***********************************************
    FluorSignal=zeros(1,F);
    for f=1:F   % for every frame
        FluorSignal(f)=mean(mean(mov(f).cdata(Meshxy_Circle(:,2),Meshxy_Circle(:,1))));
    end
    FS=[FS;FluorSignal];
    disp(['Signal: ',num2str(n),'/',num2str(NS)])
end
