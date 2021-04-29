function price = CrankNicolsonFD_BlackScholes_func(S_0, K)
% Description: Crank-Nicolson PDE Finite Difference method to price European Option in Black Scholes Model
% Author: Justin Kirkby
Smax = 3120; %max option price
Smin = 0; %min option price

N = 1600; %time steps
M = 312; %no of discrete steps in option prices

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

% z = zeros(M-1,1);
% Solve systems (backward in time)
for j = N:-1:1
%     z(1) = -a(2) * vals(1,j);
    vals(2:M,j) = U \ (L \ (M2*(vals(2:M,j+1))));
end

price = interp1(vS, vals(:,1), S_0);

figure(3);
plot(S(1:200),vals(1:200,1), 'r-', S(1:200),vals(1:200,round(N/2)), 'g-', S(1:200),vals(1:200,N+1), 'b-');
legend({'Option price at time 0', 'Option price midway through to maturity','Option price at maturity T'}, 'location', 'northwest');
xlabel('Stock price');
ylabel('Call option price');
title('European Call option using Crank-Nicholson method');

%3D plot of value of option
figure(4);
mesh(tau, S(1:310), vals(1:310,:));
zlim([0 1800]);
title('European Call Option using the Crank-Nicholson Method');
xlabel('{\tau}');
ylabel('Stock Price');
zlabel('Option Value');
end