%% Script to define Default Settings
% If Yes: Modify each part
% If No: Ask if load initial default Settings
% If Cancel: Nothing Done
Modify = questdlg('Set Default Settings for Processing?');
clc;
if strcmp('Yes',Modify)
    %% Set Directories
    % Read Experiments ***************************************
    button = questdlg('Set Experiments Directory?');
    if strcmp('Yes',button)
        DefaultPath = uigetdir(pwd);
        disp('Experiments Folder is Updated')
    else
        disp('Experiments Folder is Already Set')
    end
    % Ask if Save Processing Log Activity 
    % @ Features Tables
    % @ Resume Tables
    %% Set Biexponential Parameters and AR process order and length
    % Setup for Auto Regressive Process Estimation
    fsbiexp=10;        % Just an arbitrary sampling frequency
    Load_Default_Values_SP;
    % L=30;              % seconds  of the Fluorophore Response
    % taus_0= [.75,2,1]; % starting values
    t_r=linspace(0,L,L*fsbiexp);
    r_bi = taus_0(3)*(exp(-t_r/taus_0(2)) - exp(-t_r/taus_0(1)));
    % display biexponential function
    biexp=figure('numbertitle','off');
    biexp.Name='Biexponential Function';
    plot(t_r,r_bi,'b','LineWidth',3);
    axis tight; grid on;
    ylabel('r(t)')
    xlabel('t [s]')
    title('$r(t)=g\cdot (e^{\frac{-t}{\tau_f}}-e^{\frac{-t}{\tau_r}})$','Interpreter','latex')
    button = questdlg('Set Biexponential Parameters?');
    if strcmp('Yes',button)
        % Display deafult Values in Boxes
        % Update data plot of the biexponential
        % Accept
        disp('Default Biexpoenetial Parameters are Updated')
    else
        disp('Default Biexpoenetial Parameters are Already Set')
    end
    % p=3;               % AR(p) initial AR order        
    close(biexp);
    %% Set Denoising Setup
    %% Set Sparse Deconvolution Sttings
    %% Set Raster-Getting Configuration
elseif strcmp('No',Modify)
    button = questdlg('Set Initial Default Settings?');
    if strcmp('Yes',button)
        %% CLEAR: if not setting but setting default
        DefaultPath=pwd;
        L=30;              % seconds  of the Fluorophore Response
        taus_0= [.75,2,1]; % starting values
        p=3;               % AR(p) initial AR order        
        disp('Initial Default Settings Set')
    else
        disp('Nothing Done')
    end
    
else
    disp('Nothing has been Modified')
end
%% Save in File .mat