u0 = 10;
a = -5;
exact = @(t) u0 * exp(a * t); % exact solution

% set initial conditions and step sizes
dtm = -2 / a;
% dtq = dtm .* (2 .^ -(2 : 4));
dtq = dtm .* (2 .^ -(0 : 0.5 : 10));
rmserror = zeros(1, length(dtq));

% run through forward each step size in dtq
for k = 1 : length(dtq)
    % initialize numerical scheme for each case    
    dt = dtq(k);
    tf = 4 / dt;
    tq = dt * (0 : tf);
    u = u0;
    num = 0; % numerator of RMSE
    den = u0 ^ 2; % denominator of RMSE

    % run RK4   
    for i = 2 : length(tq)
        k1 = a * u;
        k2 = a * (u + 0.5 * dt * k1);
        k3 = a * (u + 0.5 * dt * k2);
        k4 = a * (u + dt * k3);
        u = u + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4);
        ut = exact(tq(i));
        % update RMSE        
        num = num + (u - ut) ^ 2;
        den = den + ut ^ 2;
    end
    rmserror(k) = sqrt(num / den);
end

% slope of loglog plot (error convergence rate)
m = sum((log(rmserror) - mean(log(rmserror))) .* (log(dtq) -...
    mean(log(dtq)))) / sum((log(dtq) - mean(log(dtq))) .^ 2);
disp(['Convergence rate for error: ', num2str(m), '.'])

% plot results
figure(1)
loglog(dtq, rmserror, 'o-')
% title('Problem 1(c) results (as given)', 'fontsize', 16)
title('Problem 1(c) results (with more step sizes)',...
    'fontsize', 16)
legend(['Slope: ', num2str(m)])
xlabel('$\log \Delta t$', 'interpreter', 'latex',...
    'fontsize', 16)
ylabel(['$\log\left(\frac{\sum_{i}{\left(u_{i}-u\left(t_{i}',...
    '\right)\right)^{2}}}{\sum_{i}{u\left(t_{i}\right)^2}}',...
    '\right)^{1/2}$'], 'interpreter', 'latex', 'fontsize', 16)