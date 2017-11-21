% Cumulative log normal distribution: 
% logcdf(x) = cdf(log(mu*x/sigma)), 
% where cdf is the cumulative distribution of the standard normal
% distribution, with mean=0 and variance=1

function val=logcdf(x)
    val = cdf(log(5*x)/0.15);
end