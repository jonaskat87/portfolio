u0 = 10;
a = -5;
exact = @(t) u0 * exp(a * t); % exact solution

% set initial conditions and step sizes
dtm = -2 / a;
dt = dtm / 2;
tf = 4 / dt;
tq = dt * (0 : tf);
uq = zeros(1, length(tq));
uq(1) = u0;

% run forward Euler time stepping
factor = 1 + dt * a;
for i = 2 : length(tq)
    uq(i) = factor * uq(i - 1);
end

% plot exact vs. numerical results
figure(1)
plot(tq, exact(tq), 'o-', tq, uq, 'o-')
title('Problem 1(a) solutions', 'fontsize', 16)
legend('Exact solution', ['Numerical solution ',...
    '($\Delta t=\Delta t_{m} / 2)$'], 'interpreter', 'latex')
xlabel('t', 'fontsize', 16)
ylabel('u(t)', 'fontsize', 16)

% plot errors (both in regular and semilog plots)
figure(2)
tiledlayout(2, 1)
nexttile
plot(tq, abs(uq - exact(tq)))
title('Problem 1(a) errors', 'fontsize', 16)
legend('Absolute error', 'interpreter', 'latex')
xlabel('t', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\left|u(t)-u_{h}(t)\right|$', 'interpreter',...
    'latex', 'fontsize', 16)

nexttile
semilogy(tq, abs(uq - exact(tq)), 'o-')
legend('Logarithm of the absolute error', 'interpreter', 'latex')
xlabel('t', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\log\left|u(t)-u_{h}(t)\right|$', 'interpreter',...
    'latex', 'fontsize', 16)

uq = uq(2 : end);
tq = tq(2 : end);
m = sum((log(abs(uq - exact(tq))) - mean(log(abs(uq - exact(tq))))) .* (tq -...
    mean(tq))) / sum((tq - mean(tq)) .^ 2);
disp(['Convergence rate for error: ', num2str(m), '.'])