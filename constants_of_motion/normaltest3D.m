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
for i = 1 : d + 1
    centroids(:, i) = sum(X(:, setdiff(1 : d + 1, i)), 2) / d;
    N(:, i) = N(:, i) / norm(N(:, i));
end
vecplus = centroids + N;
vecminus = centroids - N;

plot3(Xp(1, :), Xp(2, :), Xp(3, :))
hold on 
scatter3(centroids(1, :), centroids(2, :), centroids(3, :))
for i = 1 : d + 1
   plot3([vecminus(1, i) vecplus(1, i)],...
       [vecminus(2, i) vecplus(2, i)],...
       [vecminus(3, i) vecplus(3, i)]) 
end
xlim([0 1])
ylim([0 1])
zlim([0 1])
hold off