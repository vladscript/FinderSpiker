% Save Figure

function SaveFigure(figure)
snam='foo'; % The name of your style file (NO extension)
s=hgexport('readstyle',snam);

%apply style sheet info
fnam='myfig.jpeg'; % your file name
s.Format = 'jpeg'; %I needed this to make it work but maybe you wont.
hgexport(gcf,fnam,s);
[file path] = uiputfile({'*.ai','Adobe Illustrator (*.ai)';...
    '*.eps','EPS file (*.eps)';'*.jpg','JPEG image (*.jpg)';...
    '*.bmp','Bitmap file (*.bmp)';'*.fig','MATLAB figure (*.fig)';...
    '*.png','Portable Network Graphics file (*.png)'},...
    'Save file name');
if path
    saveas(figure,[path file])
end