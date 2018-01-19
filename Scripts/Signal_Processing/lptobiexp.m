% Find the taus of the biexponential funtion 
% that fit a linear predictor
% computed from linear prediction
% Input
%  tau_0 & tau_0: vector of starting taus
%  r_data: vector from the response function computed from lp
%  t_data: vector of time
% Output
%  taus, paramteres that fit the biexponential the r_data

% r = e^(-t/tau2) - e^(-t/tau1)

function [taus,model]= lptobiexp(tau_0,t_data,r_data)

% Call fminsearch
start_taus = [tau_0(1),tau_0(2),tau_0(3)];
model = @biexp;
taus = fminsearch(model, start_taus);

    function [sse, FittedCurve] = biexp(params)
        tau1 = params(1);
        tau2 = params(2);
        A = params(3);
        FittedCurve  = A*( exp(-t_data/tau2) - exp(-t_data/tau1) );
%         FittedCurve = A .* exp(-lambda * xdata);
        ErrorVector = FittedCurve - r_data;
        sse = sum(ErrorVector .^ 2); % ? mean?
    end
end