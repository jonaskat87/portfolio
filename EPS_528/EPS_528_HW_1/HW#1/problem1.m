xq = 0 : 0.001 : 100;
top = 2 * xq .^ 3;
kq = top ./ (xq .^ 2 - 1);
rq = top ./ ((xq .^ 2 + 1) .^ 2);

figure(1)
plot(kq, rq)
xticks(0 : 5 : 50)
title('Bifurcation diagram of the insect outbreak model',...
    'FontSize', 16)
xlim([0 50])    
ylim([0 0.8])
xlabel('k*(x)', 'FontSize', 16)
ylabel('r*(x)', 'FontSize', 16)
text(15, 0.7, 'outbreak', 'FontSize', 16, 'Color', 'r')
text(25, 0.35, 'bistable', 'FontSize', 16, 'Color', 'm')
text(10, 0.1, 'controlled', 'FontSize', 16, 'Color', 'b')

k = 0 : 0.02 : 15;
xq = (0 : 0.02 : 15).';
M = length(k);
N = length(xq);
kq = ones(N, 1) * k;
rq = (kq ./ (kq - xq)) .* (xq ./ (1 + (xq .^ 2)));
xq = xq * ones(1, M);

figure(2)
h1 = surf(kq, rq, xq);
set(h1,'LineStyle','none')
title('Full bifurcation diagram of the insect outbreak model',...
    'FontSize', 14)
xlim([0.01 15])
xticks(0 : 3 : 15)
ylim([0.01 0.8])
yticks(0 : 0.2 : 0.8)
zlim([0 15])
zticks(0 : 3 : 15)
xlabel('k', 'FontSize', 14)
ylabel('r', 'FontSize', 14)
zlabel('x*', 'FontSize', 14)
hold on
h2 = surf(ones(N, 1) * k, (0 : (10 / (N - 1)) : 10).' *...
    ones(1, M), zeros(N, M));
set(h2,'LineStyle','none')
hold off