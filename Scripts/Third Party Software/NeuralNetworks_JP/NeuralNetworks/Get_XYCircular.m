%% Get XY Circular
function XYCirc=Get_XYCircular(N,radio)
    if (nargin==1)
        radio=1;
    end
    
    if (N==1)
        XYCirc=[0 0];
    else
        XYCirc=[cos(2*pi*[1:N]'/N) sin(2*pi*[1:N]'/N)].*radio;
    end
end