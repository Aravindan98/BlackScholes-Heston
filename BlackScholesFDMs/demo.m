close all;
clear;
stock_data = csvread('ASIAN_PAINT_Stock_Prices.csv'); %stock level data
% 
% niftycsv = readtable('NIFTY 50_Data.csv');
% dates = datetime(niftycsv.Date,'InputFormat','dd MM yyyy');
% 
% ohlc = zeros(74,5);
% ohlc(:,1)=str2double(niftycsv.Open);
% ohlc(:,2)=str2double(niftycsv.High);
% ohlc(:,3)=str2double(niftycsv.Low);
% ohlc(:,4)=str2double(niftycsv.Close);
% ohlc(2:74,5)=diff(log(ohlc(:,4)));
S0 = 1500;
size_data = size(stock_data); %size of data
strike_prices = [1500;1520;1540;1560;1580;1600;1620;1640;1700;1760;1800];
close_price = stock_data(:,5); %select close price of last one month(22 days)

%find historical volatility by first calculating log of daily returns
vol_table = zeros(size_data(1), 2);
vol_table(:,1) = close_price;
for j=2:size_data(1)
    vol_table(j,2) = log(vol_table(j,1)/vol_table(j-1,1));
end

mean_return = mean(vol_table(:,2));
std_dev = sqrt(1/((size_data(1)-2))*sum((vol_table(2:size_data(1),2)-mean_return).^2));
sigma = std_dev*sqrt(252);
time = 26/252;
disp([std_dev,sigma, time]);
risk_free_rate = 0.0575; 

K=strike_prices(2);


Smax = 3000; %max option price
Smin = 0; %min option price

N = 1000; %time steps
M = 1000; %no of discrete steps in option prices

[sigma, r, T] = calculate_parameters();

dt = T/N; %time step
dS = (Smax - Smin)/M; %pricestep

% Description: Crank-Nicolson PDE Finite Difference method to price European Option in Black Scholes Model
% Author: Justin Kirkby
Smax = 3000; %max option price
Smin = 0; %min option price

N = 100; %time steps
M = 50; %no of discrete steps in option prices

[sigma, r, T] = calculate_parameters();

dt = T/N; %time step
dS = (Smax - Smin)/M; %pricestep

tau = dt*(0:N);

vals = zeros(M+1, N+1);
vS = linspace(Smin, Smax, M+1)';
S=vS;
vI = vS / dS;
vJ = 0:N;

% Boundary Conditions
% call option
% if call==1
vals(:, N+1) = max(vS - K, 0);  
vals(1, :) = 0;
vals(M+1, :) = Smax - K*exp(-r*dt*(N - vJ)); 
% else % put option
%     vals(:, N+1) = max(K - vS, 0);  
%     vals(1, :) = K*exp(-r*dt*(N - vJ));
%     vals(M+1, :) = 0; 
% end

% Tridiagonal Coefficients
a = 0.25 * dt * (sigma^2 *(vI.^2) - r*vI);
b = -dt*0.5*(sigma^2*(vI.^2) + r);
c = 0.25*dt*(sigma^2*(vI.^2) + r*vI);

M1 = -diag(a(3:M), -1) + diag(1 - b(2:M)) - diag(c(2:M-1),1);
[L,U] = lu(M1);
M2 = diag(a(3:M), -1) + diag(1 + b(2:M)) + diag(c(2:M-1),1);
% Solve systems (backward in time)
for j = N:-1:1
%     vals(2:M,j) = U \ (L \ (M2*vals(2:M,j+1)));
    vals(2:M,j) = pinv(M1)*(M2*(vals(2:M,j+1)+z));
end

price = interp1(vS, vals(:,1), S0+1500);

figure(3);
plot(S,vals(:,1), 'r-', S,vals(:,round(N/2)), 'g-', S,vals(:,N+1), 'b-');
legend({'Option price at time 0', 'Option price midway through to maturity','Option price at maturity T'}, 'location', 'northwest');
xlabel('Stock price');
ylabel('Call option price');
title('European Call option using Crank-Nicholson method');

%3D plot of value of option
figure(4);
mesh(tau, S, vals);
zlim([0 1700]);
title('European Call Option using the Crank-Nicholson Method');
xlabel('{\tau}');
ylabel('Stock Price');
zlabel('Option Value');