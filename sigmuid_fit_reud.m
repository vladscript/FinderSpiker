function [taus,model]= sigmuid_fit_reud(Coeff,r_data)
t_data=1:numel(r_data);
% Call fminsearch
start_taus = [Coeff(1),Coeff(2)];
model = @sigmuid;
taus = fminsearch(model, start_taus);

    function [sse, FittedCurve] = sigmuid(params)
        a = params(1);
        b = params(2);
        FittedCurve  = a./ ( 1+ exp(-b*t_data));
        ErrorVector = FittedCurve - r_data;
        sse = sum(ErrorVector .^ 2); % ? mean?
    end
end