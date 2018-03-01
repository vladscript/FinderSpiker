% Funtion to Read Images from another lambda
% and find cell that colocates with Calcium Imaging Experiments
% It reads any image and display with the population coordinates
% It allows to the user determine if it colocates (automatization...)
% Save Colocated Coordinates in the .mat File of the Experiment
% Green Coordinates: Colocated Neurons
% Red Coordinates: No Colocated Neurons
% Input
%  Experiment:  Experiment ID
%   XY:         Population Coordinates
%   r:          Radious Lengtth
%               Get File using Dialog
% 
% Ouput
%               Display Image and Coordinates
% XY_merged     Colocated Coordinates
function [XY_merged,ColocateIndx]=get_merged_coordinates(Experiment,XY,r)
%% Setup
N=length(XY);
XY_merged=[];
X_col=[];
Y_col=[];
% ColocateColor=['r','g'];
ColocateIndx=zeros(N,1);

%% Open Dialogue Box to Read Image
DefaultPath='C:\Users\Vladimir\Documents\Doctorado\Experimentos';
        if exist(DefaultPath,'dir')==0
            DefaultPath=pwd; % Current Diretory of MATLAB
        end
[FileNameImage,PathNameImage] = uigetfile('*.bmp',[' Pick the Analysis File ',Experiment],...
            'MultiSelect', 'off',DefaultPath);
[A,map]=imread([PathNameImage,FileNameImage]);  % Images
B=ind2gray(A,map);                              % Gray Scale Image
[H,W]=size(B);                                  % Size of the Image
%% PLot Figure
getcoord=figure('numbertitle','off',...
    'name',['Colocalization @',Experiment]);

%     'keypressfcn',@exit_function,...
%     'ButtonDownFcn',@selec_ROI);
%% Show Image
imshow(B);
imcontrast;
hold on;
% Plot Coordinates & Get ROI pixels
Meshxy_Circle=[]; 
aux1=1;
for n=1:N
%     rectangle('Position',[XY(n,1)-r(n)/2,XY(n,2)-r(n)/2,2*r(n),2*r(n)],...
%             'Curvature',[1 1],'LineWidth',1,'EdgeColor','r');
    rectangle('Position',[XY(n,1)-r(n),XY(n,2)-r(n),2*r(n),2*r(n)],...
        'Curvature',[1 1],'LineWidth',1,'EdgeColor','r');
    % ROI pixels
    Mx=XY(n,1)-(r):XY(n,1)+(r); % range in x of square
    My=XY(n,2)-(r):XY(n,2)+(r); % range in y of square
    for i=1:length(Mx)
        % chech if it's in image's limits Xaxis
        if Mx(i)>0 && Mx(i)<=W 
            for j=1: length(My)
                % chech if it's in image's limits Yaxis
                if My(j)>0 && My(j)<=H 
                    % check if it's in circle
                    if (Mx(i)-XY(n,1))^2+(My(j)-XY(n,2))^2<=r(n)^2
                        Meshxy_Circle=[Meshxy_Circle;Mx(i),My(j)];
                        ROIcounter(n)=aux1;
                        aux1=aux1+1;
                    end
                end
            end
        end
    end
end
% WAIT TO SEE IMAGE
% pause;
%% Get Points 
% IT GETS POINTS ONLY IN THE CIRCLES OF THE COORDINATES
Meshxy_Circle=round(Meshxy_Circle);
[X_col,Y_col]=getpts(getcoord);
X_col=round(X_col);
Y_col=round(Y_col);
% CHeck POint are in ROIs
N_col=length(X_col);
AcceptedIndx=[];
for i=1:N_col
    [~,location]=ismember([X_col(i),Y_col(i)],Meshxy_Circle,'rows');
    % if ismember([X_col(i),Y_col(i)],Meshxy_Circle,'rows')
    % if ~isempty(location)
    if ~location==0
        AuxInd=find(ROIcounter>location);
        ColocateIndx(AuxInd(1))=1;
        XY_merged=[XY_merged;XY(AuxInd(1),:)];
        AcceptedIndx=[AcceptedIndx,i];
    end
end

disp('Rejected: ')
disp(N_col-length(AcceptedIndx))
if isempty(AcceptedIndx)
    disp('***********NO CO-LOCALIZATION***************')
    hold off;
    imshow(B);
else
    hold off;
    imshow(B);
    imcontrast;
    hold on;
    % Plot Coordinates & Get ROI pixels
    for n=1:length(AcceptedIndx)
        rectangle('Position',[XY_merged(n,1)-r(n),XY_merged(n,2)-r(n),2*r(n),2*r(n)],...
        'Curvature',[1 1],'LineWidth',1,'EdgeColor','g');
    end
    % Save IMAGE Manually
end



% X_col=X_col(AcceptedIndx);
% Y_col=Y_col(AcceptedIndx);
% % Save Colocated Coordinates
% if ~isempty(Y_col)
%     N_col=length(X_col);
%     for i=1:N_col
%         for j=1:N
%             dis(j)=abs(X_col(i)-XY(j,1))+abs(Y_col(i)-XY(j,2));
%         end
%         [~,locind]=min(dis);
%         ColocateIndx(locind)=1;
%         XY_merged=[XY_merged;XY(locind,:)];
%     end
% else
%     disp('NO CO-LOCALIZATION')
% end

%% Save Colocated Coordinates
checkname=1;
while checkname==1
    DefaultPath='C:\Users\Vladimir\Documents\Doctorado\Software\GetTransitum\Calcium Imaging Signal Processing\FinderSpiker\Processed Data';
    if exist(DefaultPath,'dir')==0
        DefaultPath=pwd; % Current Diretory of MATLAB
    end
    [FileName,PathName] = uigetfile('*.mat',[' Pick the Analysis File ',Experiment],...
        'MultiSelect', 'off',DefaultPath);
    dotindex=find(FileName=='.');
    if strcmp(FileName(1:dotindex-1),Experiment(2:end))
        checkname=0;
        % SAVE DATA
        save([PathName,FileName],'XY_merged','ColocateIndx','-append');
        disp([Experiment,'   -> UPDATED (Merged Coordinates)'])
    else
        disp('Not the same Experiment!')
        disp('Try again!')
    end
end    

%% Nested 
%     function exit_function(getcoord,~,~)
%        key=get(getcoord,'CurrentKey');
%        if strcmp(key,'e') % Exit
%             close(getcoord);
%        end
%     end
 
%     function selec_ROI(~,~)
%         [X_col,Y_col]=getpts;
%         analyze_colocalization();
%     end

%     function analyze_colocalization()
%         N_col=length(X_col);
%         for i=1:N_col
%             for j=1:N
%                 dis(j)=abs(X_col(i)-XY(j,1))+abs(Y_col(i)-XY(j,2));
%             end
%             [~,locind]=min(dis);
%             ColocateIndx(locind)=1;
%             XY_merged=[XY_merged;XY(locind,:)];
%         end
%     end
end

%% END OF THE WORLD