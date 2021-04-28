function call_price = Black_Scholes_Script(S0, E)

% Reference[1] for implicit method is
% https://www.math.uaic.ro/~annalsmath/pdf-uri%20anale/F1(2010)/Mosneagu.pdf(Numerical
% Approximation of Black Scholes Equation)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 1600; %time steps
N = 312; %no of discrete steps in option prices
% Ensure that M > N, else this is not implicitly uniquely solvable.
Smax = 3120; %max option price
Smin = 0; %min option price
[sigma, r, T] = calculate_parameters();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dT = T/M; %time step
dO = (Smax - Smin)/N; %pricestep

%initialise matrices
V = zeros(N+1, M+1); % option prices 
S = Smin + dO*(0:N); % list of option price steps
tau = dT*(0:M); %list of time steps

% Initial condition is basically final condition,
% so we use tau = T - t instead of t and condition is 
% V(S, tau=0) = max(S-E,0)

V(1:N+1, 1) = max(S-E,0);

%Boundary conditions:
% 1. V(0,t) = 6
% 2. V(S,t) = S - E*exp(-r*tau) as S -> inf or when S = Smax
% Reference[2] for BCs and IC:
% http://math.yorku.ca/~dhackman/BlackScholes7.pdf(Solving the Black Scholes Equation using a
% Finite Difference Method)

V(1,1:M+1) = 0.0; %BC 1
V(N+1,1:M+1) = Smax - E*exp(-r*tau); %BC 2

% Note that this vectorised implementation is not part of the report. Please
% refer Dura, Gina "Numerical Approximation of Black Scholes Equation" for
% the vectorised implementation under explicit method. 
% Reference[1] is
% https://www.math.uaic.ro/~annalsmath/pdf-uri%20anale/F1(2010)/Mosneagu.pdf
% Refer sec 2.2 in Reference[1] for terms alpha beta and finite difference formulation
alpha = (sigma^2)*dT;
beta = r*dT;

%Everything is according to sec 2.2 of Reference[1]- explicit method
%formulation
%Initialise A matrix and other variables
A = zeros(N-1,N-1);



%Fill A matrix 
for i = 1:N-1
    if i ~= 1
        A(i-1,i) = -0.5*(alpha*((i-1)^2) + beta*(i-1)); % un
    end
    A(i,i) = 1 + alpha*(i^2) + beta; % dn 
    if i+1 ~= N
        A(i+1,i) = 0.5*(beta*(i+1) - alpha*((i+1)^2)); % ln 
    end
end

l0 = 0.5*(beta-alpha);
un = -0.5*(alpha*((N-1)^2) + beta*(N-1));
Z = zeros(N-1,1);

for j = 2:M+1
    Z(1) =  l0*V(1,j);
    Z(N-1) = un*V(N+1,j);
    B = V(2:N,j-1)-Z;
    V(2:N,j) = A\B;
    %t1 = linsolve(A,B);
    % The system to solve is A*V(2:N,j)=B to find next values of V(2:N,j)
end

%Let's figure out option price to return
index = round((S0 - Smin)*N/(Smax - Smin))+ 1;
call_price = V(index, M+1);

%Figure of value of option, V(S, tau) as a function of S at three diff
%times: tau = 0(T= t), tau= T/2(t = T/2) and tau = T(t = 0)
figure(1)
plot(S,V(:,1), 'r-', S,V(:,round(M/2)), 'g-', S,V(:,M+1), 'b-')
legend({'Option price at maturity T', 'Option price midway through to maturity', 'Option price at time 0'}, 'location', 'northwest')
xlabel('Stock price')
ylabel('Call option price')
title('European Call option using implicit method')

%3D plot of value of option
figure(2)
mesh(tau, S, V)
zlim([0 1700])
title('European Call Option using the Implicit Method')
xlabel('{\tau}')
ylabel('Stock Price')
zlabel('Option Value')