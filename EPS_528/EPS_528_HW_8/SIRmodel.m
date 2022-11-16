function [ ] = SIRmodel(eps, mu, gamma, b1, b2, b3, Omega1, ...
    I0, R0, T)

S0 = 1 - R0 - I0;
beta = @(t)  16 * mu * (1 + b1 * (eps ^ 2) + b2 * (eps ^ 4) +...
    (eps ^ 5) * b3 * sin(mu * eps * Omega1 * t));
SIR = @(t, y) [mu * (1 - y(1)) - beta(t) * (y(2) ^ 2) * y(1);...
    beta(t) * (y(2) ^ 2) * y(1) - (gamma + mu) * y(2);
    gamma * y(2) - mu * y(3)];

tol = 1e-10;
opts = odeset('RelTol', tol, 'AbsTol', tol);
[t, y] = ode45(SIR, [0 T], [S0; I0; R0], opts);

figure
plot(t, y)
legend

figure
scatter3(y(:, 1), y(:, 2), y(:, 3))
xlabel('S')
ylabel('I')
zlabel('R')