% Gets Upper Directory
function FileDirSave=getDir2Save()
FileDirSave=pwd;
slashes=find(FileDirSave=='\');
FileDirSave=FileDirSave(1:slashes(end));