function [X,fs] = readdeltaFzero(filename)

X=table2array(readtable(filename,'ReadVariableNames',false,'HeaderLines',4))';
delimiter = ',';
fileID = fopen(filename,'r');
fsamples = textscan(fileID, '%f',1, 'Delimiter', delimiter, 'HeaderLines', 0, 'ReturnOnError', false);
fclose(fileID);
fs=fsamples{1};