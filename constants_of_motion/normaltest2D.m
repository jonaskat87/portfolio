d = 2; % do not change
X = rand(d, d + 1);
A = [1 2];
B = [2 3];
C = [3 1];
Xp = X(:, [A B C]);
Einv = inv([transpose(X) ones(d + 1, 1)]);
N = Einv(1 : end - 1, :);
centroids = zeros(size(X));
for i = 1 : d + 1
    centroids(:, i) = sum(X(:, setdiff(1 : d + 1, i)), 2) / d;
    N(:, i) = N(:, i) / norm(N(:, i));
end
vecplus = centroids + N;
vecminus = centroids - N;

plot(Xp(1, :), Xp(2, :))
hold on 
scatter(centroids(1, :), centroids(2, :))
for i = 1 : d + 1
   plot([vecminus(1, i) vecplus(1, i)],...
       [vecminus(2, i) vecplus(2, i)]) 
end
xlim([0 1])
ylim([0 1])
hold off