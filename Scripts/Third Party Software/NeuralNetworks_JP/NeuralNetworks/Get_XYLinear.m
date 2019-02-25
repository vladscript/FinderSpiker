%% Get XY Linear
function XYLine=Get_XYLinear(N)
    XYLine=[repmat(-1.1,N,1) ((1:N)*2/N)'-1];
end