%This is the main script
%We are looking at the simulated option prices and corresponding premium of
%Asian Paints Stock Option
%Following is the historical option price data ref: 
%https://fnoanalysis.com/oi/option_chain_hist.php?symbol=ASIANPAINT&cmb_cnd_symbol=ASIANPAINT&CMB_EXPIRY_DT=2020-06-25&CMB_CND_DT=2020-05-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We are simulating call option price for different strike prices and
%comparing it to the premium or running charge of buying the option
S0 = 1500; %stock price on 20th May
S_maturity = 1688.85; %stock price on 25th June
strike_prices = [1500;1520;1540;1560;1580;1600;1620;1640;1700;1760;1800];
payoff = S_maturity - strike_prices;
a=-0.030; %Lower bound of the uniform jump-amplitude mark density 
b=0.028; %Upper bound of the uniform jump-amplitude mark density
lambda=40; %lambda  - Annualized jump rate 
m = 1000; % Number of Monte Carlo Simulations
call_implicit = zeros(11,1);
call_exact = zeros(11,1);
call_explicit =zeros(11,1);
call_CN = zeros(11,1);
abs_error_exp = zeros(11,1);
abs_error_imp = zeros(11,1);
abs_error_CN = zeros(11,1);


premium = [74.85;65.00;60.00;51.40;57.20;36.65;35;37.40;19.80;15.00;11.10];

for i = 1:11
    %call_calc(i) = CrankNicolsonFD_BlackScholes_func(S0,strike_prices(i));
    call_implicit(i) = Black_Scholes_Script(S0, strike_prices(i));
    call_exact(i) = exact_call_price(S0,strike_prices(i));
    call_explicit(i) = Explicit_BS(S0,strike_prices(i));
    call_CN(i) = CrankNicolsonFD_BlackScholes_func(S0,strike_prices(i));
    abs_error_imp(i) = abs(call_exact(i) - call_implicit(i))/call_exact(i);
    abs_error_CN(i) = abs(call_exact(i) - call_CN(i))/call_exact(i);
    abs_error_exp(i) = abs(call_exact(i) - call_explicit(i))/call_exact(i);
end

call_analytical=call_exact;
analysis_table = table(strike_prices, call_analytical, call_explicit, call_implicit, call_CN);
disp(analysis_table);
analysis_table.Variables = round(analysis_table.Variables,5);
writetable(analysis_table,'price.csv','Delimiter',',');


error_table = table(strike_prices,abs_error_exp, abs_error_imp, abs_error_CN);
disp(error_table);
error_table.Variables = round(error_table.Variables,5);
error_table.Variables = round(error_table.Variables,5);
writetable(error_table,'error.csv','Delimiter',',');

plot(1:11,abs_error_exp,'r',1:11,abs_error_imp,'g',1:11,abs_error_CN,'b');
legend({'explicit','implicit','crank-nicholson'});
title('Absolute Relative Error');
