function [ data ] = percolation3D(N, dp)
tic
format long

pq = 0 : dp : 1;
Mq = zeros(size(pq));
for i = 1 : length(pq)
   % generate N 2x2x2 sample cubes
   samples = rand(2, 2, 2, N) > 1 - pq(i);
   % sum over top-to-bottom paths
   % (sum will be 2 if path exists in that column)
   paths = squeeze(sum(samples, 3));
   % check for paths
   checkpaths = reshape(paths, [4, N]) == 2;
   % if contains path, sum of column will be greater than
   % 0. use this to count the number of permeable cubes 
   Mq(i) = sum(sum(checkpaths) > 0);
end
Mq = Mq / N;

R2 = @(p) 4 * p .^ 2 - 6 * p .^ 4 + 4 * p .^ 6 - p .^ 8;
R2zero = @(p) R2(p) - p;
% 0.3 is close to the p^{*}, so fzero should converge to 
% p^{*} if we let 0.3 be our initial guess for fzero
pc = fzero(R2zero, 0.3);

% root-mean-square error between sampling 
% and theoretical prediction
RMSE = sqrt((R2(pq) - Mq) * ((R2(pq) - Mq).') / length(pq));

figure
hold on
plot(pq, pq, ':', 'Color', 'k', 'HandleVisibility', 'off') 
plot(pq, Mq, '.', 'Color', '#77AC30')
plot(pq, R2(pq), 'Color', '#4DBEEE')
plot(0, 0, 'r*')
plot(pc, R2(pc), 'm*') 
plot(1, 1, 'b*')
hold off
grid on
legend({['Numerical (RMSE: ', num2str(RMSE), ')'], ...
    'Theory', '$\left(0,0\right)$',...
    ['$\left(p^{*}, R_{2}\left(p^{*}\right)\right)\approx',...
    '\left(', num2str(pc), ',', num2str(pc), '\right)$'],...
    '$\left(1,1\right)$', }, 'Location', 'southeast',...
    'interpreter', 'latex', 'FontSize', 12)
title(['3-D site percolation simulation with $N=$',...
    num2str(N), ' and $\Delta p=$', num2str(dp)],...
    'interpreter', 'latex', 'FontSize', 16)
xlabel('$p$ (probability of permeability)',...
    'interpreter', 'latex', 'FontSize', 14)
ylabel(['$R_{2}(p)\approx M/N$ (renormalization',...
    ' group transformation)'], 'interpreter', 'latex',...
    'FontSize', 14)

dR2 = 8 * pc - 24 * pc .^ 3 + 24 * pc .^ 5 - 8 * pc .^ 7;
v = log(2) / log(dR2);

disp('The left value is the critical occupation probability.')
disp('The right value is the critical exponent.')
data = [pc v];

toc