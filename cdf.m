% Cumulative distribution of the standard normal
% distribution, with mean=0 and variance=1
function val=cdf(x)
    val = 0.5*(1+erf(x/sqrt(2.0)));
end