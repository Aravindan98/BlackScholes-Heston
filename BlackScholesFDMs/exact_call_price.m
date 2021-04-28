function call = exact_call_price(S0, strike)
% Parameters
% Need to find an actual options and parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%r = 0.2 ; %risk free rate
%sigma = 0.25; %volatility
%T = 0.5; time period
%Above are only examples. Not to be used.
S = S0;
E = strike;
[sigma, r, T] = calculate_parameters();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1 = (log(S/E) + (r +0.5*sigma*sigma)*(T))/(sigma*sqrt(T));
d2 = (log(S/E) + (r -0.5*sigma*sigma)*(T))/(sigma*sqrt(T));
call = S*normcdf(d1) - E*exp(-r*T)*normcdf(d2);