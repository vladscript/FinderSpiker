function Modo=SelectMode()
%% Initialize Figure
global Modo;
Modo='none';
fprintf('\n>>Select FinderSpiker mode:')
IntroFig=figure('DockCOntrols','off','MenuBar','none', 'NumberTitle',...
    'off','ToolBar','none','Units','normalized','Position',...
    [0.3 0.3 0.4 0.25],'Name','FinderSpiker.CalciumImaging','Visible','off');
imshow('Figures\Logo_FinderSpiker.png')
ButtonM=uicontrol('Parent',IntroFig,'Style','pushbutton','String',...
    'ImPatch Analysis','Units','normalized','Position',[0.0 0.0 0.4 0.2],...
    'Visible','on','Callback',@ActionButtonM);
ButtonS=uicontrol('Parent',IntroFig,'Style','pushbutton','String',...
    'CSV Analysis','Units','normalized','Position',[0.6 0.0 0.4 0.2],...
    'Visible','on','Callback',@ActionButtonS);
IntroFig.Visible='on';

waitfor(IntroFig);
%% Functions
function ActionButtonM(~,~)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf(' ImPatch mode\n');
Modo='Multiple';
close(IntroFig);
end
function ActionButtonS(~,~)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fprintf(' Single CSV mode\n');
Modo='Single';
close(IntroFig);
end
end