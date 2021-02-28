% Function to plot Rainclouds for multilple paired/unpaired measures
% Input
%   MMeasures:  String 'paired' onnly if paired measures
%   ColorsN:    Matrix Colors for N inputs,
%   bda_init:   Plot parameter,
%   bdam_init   Plot parameter,
%   varargin:   N vector Inputs or Cell
%  Output
%   bda_final   Plot parameter for following boxes
%   bdam_final  Plot parameter for following boxes
% 
% Example 1
% Three repeated measures: A,B,C
% figure
% >>[bda_final,bdam_final]=pairedraincloud('paired',ColorsN,0,0,A,B,C)
% 
% Example 2
% Three independent measures: A,B,C
% figure
% >>[bda_final,bdam_final]=pairedraincloud('unpaired',ColorsN,0,0,A,B,C)
% 
% Example 3:
% Two experiments: A,B, X,Y
% figure
% >>[bda_final,bdam_final]=pairedraincloud('paired',ColorsN,0,0,A,B)
% >>[bda_final,bdam_final]=pairedraincloud('paired',ColorsN,bda_final,bdam_final,X,Y)


function [bda_final,bdam_final]=pairedraincloud(MMeasures,ColorsN,bda_init,bdam_init,varargin)
%% Setup
Nc=numel(varargin);
StepDodge=0.8;
Step_bdam=0.23;
box_dodge_amount=bda_init+StepDodge;
bdam=bdam_init+Step_bdam;
% Check if is Cell input
fprintf('>Input data Type: ')
if Nc==1 && iscell(varargin(1))
    Data2plot=varargin{1};
    Nc=numel(Data2plot);
    fprintf('cell.\n')
else    
    Data2plot=varargin;
    fprintf('vectors.\n')
end
    
%% Rainclouds
for n=1:Nc
    raincloud_plot(Data2plot{n},'color',ColorsN(n,:),'box_on',1,...
        'alphaval',5,'box_dodge',5,...
        'box_dodge_amount',box_dodge_amount, 'dot_dodge_amount', box_dodge_amount,...
        'box_col_match',0,'box_dodge_amount',bdam, 'dot_dodge_amount', bdam,...
        'line_width',3,'lwr_bnd',2);
    if n<Nc
        box_dodge_amount=box_dodge_amount+StepDodge;
        bdam=bdam+Step_bdam;
    end
end
bda_final=box_dodge_amount;
bdam_final=bdam;
drawnow;
%% Line for Paired/Repeated Measures
if strcmp(MMeasures,'paired')
    fprintf('>>Repeated Measures:');
    Bar1=gca;
    h1=findobj(Bar1,'Type','Scatter');
    
    if numel(h1)==Nc
        fprintf(' %i OK',Nc);
    end
    
    for n=2:Nc
        Xdata1=h1(n-1).XData;
        Ydata1=h1(n-1).YData;
        Xdata2=h1(n).XData;
        Ydata2=h1(n).YData;
        for i=1:numel(Data2plot{n}(:,1))
            plot([Xdata1(i),Xdata2(i)],[Ydata1(i),Ydata2(i)],...
                'LineStyle','-','Color','k','LineWidth',0.1)
        end
    end
else
    disp('>>Repeated Measures');
end