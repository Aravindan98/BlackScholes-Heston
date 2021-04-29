function call_price = Explicit_BS(S_0, K)

N = 1600; % time steps
M = 312; % no of discrete steps in option prices
% Ensure that M > N, else this is not implicitly uniquely solvable.
Smax = 3120; % max option price
Smin = 0; % min option price
[sigma, r, T] = calculate_parameters();
dO = (Smax - Smin) / M; % pricestep
S = Smin + dO * (0: M); % list of option price steps

% dS = Smax / M; % readjust
dt = T / N;  % readjust
tau = dt * (0: N); % list of time steps

vals = zeros(M + 1, N + 1);
vS = linspace(Smin, Smax, M + 1)';
     vI = 0:M;
     vJ = 0:N;

     call = 1;
     % Boundary Conditions
     if call == 1  % call option
     vals(:, N+1) = max(vS - K, 0);
     vals(1, :) = 0;
     vals(M+1, :) = Smax - K*exp(-r*dt*(N - vJ));
     else % put option
     vals(:, N+1) = max(K - vS, 0);
     vals(1, :) = K*exp(-r*dt*(N - vJ));
     vals(M+1, :) = 0;
     end

     % Tridiagonal Coefficients
     a = 0.5 * dt * (sigma^2*vI - r).*vI;
     b = 1 - dt * (sigma^2*vI.^2 + r);
     c = 0.5 * dt * (sigma^2*vI + r).*vI;

     % Solve (backward)
     for j = N:-1:1
     for i=2:M
     vals(i,j) = a(i)*vals(i-1,j+1) + b(i)*vals(i,j+1) + c(i)*vals(i+1,j+1);
     end
     end

     call_price = interp1(vS, vals(:,1), S_0);


     %Let's figure out option price to return

         % Figure of value of option, V(S, tau) as a function of S at three diff
         % times: tau = 0(T = t), tau = T / 2(t = T / 2) and tau = T(t = 0)
                                        figure(1)
                                        plot(S, vals(:, 1), 'r-', S, vals(:, round(N / 2)), 'g-', S, vals(:, N + 1), 'b-')
                                        legend({'Option price at maturity T', 'Option price midway through to maturity', 'Option price at time 0'}, 'location', 'northwest')
                                        xlabel('Stock price')
                                        ylabel('Call option price')
                                        title('European Call option using the Explicit method')

                                        % 3D plot of value of option
                                        figure(2)
                                        mesh(tau, S, vals)
                                        zlim([0 1700])
                                        title('European Call Option using the Emplicit Method')
                                        xlabel('{\tau}')
                                        ylabel('Stock Price')
                                        zlabel('Option Value')

                                        end