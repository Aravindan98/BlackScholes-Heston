function price = ImplicitFD_BlackScholes_func(S_0, K)
% Description: Fully Implicit PDE Finite Difference method to price European Option in Black Scholes Model
% Author: Justin Kirkby
Smax = 3000; %max option price
Smin = 0; %min option price

N = 1600; %time steps
M = 300; %no of discrete steps in option prices

[sigma, r, T]=calculate_parameters();

dt = T/N; %time step
dS = (Smax - Smin)/M; %pricestep

tau = dt*(0:N);  % readjust

vals = zeros(M+1, N+1);
vS = linspace(Smin, Smax, M+1)';
S=vS;
vI = 0:M;
vJ = 0:N;

% Boundary Conditions
 % call option
vals(:, N+1) = max(vS - K, 0);  
vals(1, :) = 0;
vals(M+1, :) = Smax - K*exp(-r*dt*(N - vJ)); 


% Tridiagonal Coefficients
a = 0.5 * (r * dt * vI - sigma^2*dt*(vI.^2));
b = 1 + sigma^2*dt*(vI.^2) + r*dt;
c = - 0.5 * (r * dt * vI + sigma^2*dt*(vI.^2));

cvec = diag(a(3:M), -1) + diag(b(2:M)) + diag(c(2:M-1),1);
[L,U] = lu(cvec);

% Solve systems (backward in time)
z = zeros(M-1,1);
for j = N:-1:1
    %z(1) = -a(2) * vals(1,j);
    vals(2:M,j) = U \ (L \ (vals(2:M,j+1) + z));
end

price = interp1(vS, vals(:,1), S_0);
figure(1);
plot(S,vals(:,1), 'r-', S,vals(:,round(N/2)), 'g-', S,vals(:,N+1), 'b-');
legend({'Option price at time 0', 'Option price midway through to maturity','Option price at maturity T'}, 'location', 'northwest');
xlabel('Stock price');
ylabel('Call option price');
title('European Call option using Implicit method');

%3D plot of value of option
figure(2);
mesh(tau, S, vals);
zlim([0 1800]);
title('European Call Option using the Implicit Method');
xlabel('{\tau}');
ylabel('Stock Price');
zlabel('Option Value');
end