tic
d = 50;
N = 5e3;
lower = -10;
upper = 10;
vec = (upper - lower) * rand(N, d) + lower;
pointsnear = zeros(d + 1, N);
vol = zeros(1, N);
volj = zeros(d + 1, N);
for i = 1 : N
    dist = vec - vec(i, :);
    dist = sum(dist .* dist, 2);
    [~, I] = sort(dist);
    ind = I(2 : d + 2);
    pointsnear(:, i) = ind;
    X = vec(ind, :);
    Einv = inv([X ones(d + 1, 1)]);
    N = Einv(1 : end - 1, :);
    X = X(2 : end, :) - X(1, :);
    G = X * transpose(X);
    vol(i) = sqrt(det(G)) / factorial(d);
    for j = 1 : d + 1
        N(:, j) = N(:, j) / norm(N(:, j));
        Xj = vec(setdiff(ind, ind(j)), :);
        Xj = Xj(2 : end, :) - Xj(1, :);
        Gj = Xj * transpose(Xj);
        volj(j, i) = sqrt(det(Gj)) / factorial(d - 1);
    end
end
toc