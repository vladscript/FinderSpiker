function [X,fs] = readdeltaFzero(filename)

X=table2array(readtable(filename,'ReadVariableNames',false,'HeaderLines',4))';
delimiter = ',';
fileID = fopen(filename,'r');

fprintf('\n> Reading CSV file: ')
fsamples = textscan(fileID, '%f',1, 'Delimiter', delimiter, 'HeaderLines', 0, 'ReturnOnError', true);
if isempty(fsamples{1})
    fprintf('raw fluorescence data.\n');
    fclose(fileID);
    fileID = fopen(filename,'r','n','UTF-8');
    % Skip the BOM (Byte Order Mark).
    fseek(fileID, 3, 'bof');
    fsamples = textscan(fileID, '%f',1, 'Delimiter', delimiter, 'HeaderLines', 0, 'ReturnOnError', true);
else
    fprintf('normallized fluorescence data.\n')
end
fclose(fileID);
fs=fsamples{1};