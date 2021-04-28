function call = BS_Analytical(S0, strike, T)
% Parameters
% Need to find an actual options and parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%r = 0.2 ; %risk free rate
%sigma = 0.25; %volatility
%T = 0.5; time period
%Above are only examples. Not to be used.
S = S0;
E = strike;
r=4.85/100;
sigma=0.00723;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d1 = (log(S/E) + (r +0.5*sigma*sigma)*(T))/(sigma*sqrt(T));
d2 = (log(S/E) + (r -0.5*sigma*sigma)*(T))/(sigma*sqrt(T));
call = S*normcdf(d1) - E*exp(-r*T)*normcdf(d2);