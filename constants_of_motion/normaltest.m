d = 3; % do not change
X = rand(d, d + 1);
A = [2 3 4 2];
B = [1 3 4 1];
C = [1 2 4 1];
D = [1 2 3 1];
Xp = X(:, [A B C D]);
Einv = inv([transpose(X) ones(d + 1, 1)]);
N = Einv(1 : end - 1, :);
centroids = zeros(size(X));
for i = 1 : 4
    centroids(:, i) = sum(Xp(:, 4 * i - 3 : 4 * i - 1), 2) / 3;
    N(:, i) = N(:, i) / (10 * norm(N(:, i)));
end
Nnew = zeros(size(N));
for i = 1 : 4
    Nnew(:, i) = Xp(:, 4 * i - 3 : 4 * i - 1) \ ones(3, 1);
    Nnew(:, i) = Nnew(:, i) / (10 * norm(Nnew(:, i)));
end
norm(N - Nnew)
vec = centroids + Nnew;

plot3(Xp(1, :), Xp(2, :), Xp(3, :))
hold on 
scatter3(centroids(1, :), centroids(2, :), centroids(3, :))
for i = 1 : 4
   plot3([centroids(1, i) vec(1, i)],...
       [centroids(2, i) vec(2, i)],...
       [centroids(3, i) vec(3, i)]) 
end
hold off