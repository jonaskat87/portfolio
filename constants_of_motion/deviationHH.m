function [sigmas, choices] = deviationHH(k, u0, ti, tf,...
    lambda, tol, snr, sketch)
% u0 = [0.1 0.4 -0.3 -0.2];
% ti = 0;
% tf = 100;
% lambda = 1;
% tol = 1e-8;
% snr = 50; % signal-to-noise ratio
% sketch = 0; % do you want to use a low-rank matrix sketch or not?

f = @(t, x, px, y, py) [px; -x * (1 + 2 * lambda * y);...
    py; -y - lambda * (x^2 - y^2)];
func = @(t, u) f(t, u(1), u(2), u(3), u(4));
options = odeset('RelTol', tol, 'AbsTol', tol);
[~, u] = ode45(func, [ti tf], u0, options);
u = awgn(u, snr);

figure(2)
plot(u(:, 1), u(:, 3))

g = @(x, px, y, py) [x px y py ...
    x .* y x .* px x .* py px .* y px .* py y .* py ...
    x .^ 2 px .^ 2 y .^ 2 py .^ 2 ...
    y .* (x .^ 2) y .^ 3];
gsize = g(1, 1, 1, 1);

n = length(u(:, 1));
M = zeros(n, length(gsize));
cumsum = g(u(1, 1), u(1, 2), u(1, 3), u(1, 4));
for i = 1 : n - 1
    term = g(u(i + 1, 1), u(i + 1, 2), u(i + 1, 3), u(i + 1, 4));
    M(i, :) = cumsum - i * term;
    cumsum = cumsum + term;
end
if sketch
    [~, ~, V] = svdsketch(M);
    sigmas = vecnorm(M * V(:, end - k + 1 : end)) / norm(M);
    choices = V(:, end - k + 1 : end) ./...
        max(abs(V(:, end - k + 1 : end)));
else
    [V, D] = eigs(transpose(M) * M, k, 'smallestabs');
    sigmas = sqrt(diag(D)) / norm(M);
    choices = V;
end