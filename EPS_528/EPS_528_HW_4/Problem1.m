cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\EPS_528_HW_4'

sigma = 10;
b = 8 / 3;
r = 28;
x0 = 0;
y0 = 1;
z0 = 0;
t0 = 0;
tf = 500;
sol = lorenz(sigma, r, b, t0, tf, x0, y0, z0);
t = sol.x;
u = sol.y;

% create Lorenz map
maxima = NaN(1, length(t));
samplei = u(3, 1);
sampleafter = u(3, 2);
for i = 2 : length(t) - 1
    samplebefore = samplei;
    samplei = sampleafter;
    sampleafter = u(3, i + 1);
    if (samplei > samplebefore) && (samplei > sampleafter)
        maxima(i) = samplei;
    end
end
maxima = maxima(~(isnan(maxima)));

figure(1)
tiledlayout(3, 1)
nexttile
plot(t, u(1, :), 'Color', [0 0.4470 0.7410])
title(['Lorenz equations solution from t=', num2str(t0),...
    ' to t=', num2str(tf)], 'fontsize', 14)
xlabel('t', 'fontsize', 12)
ylabel('x(t)', 'fontsize', 12)

nexttile
plot(t, u(2, :), 'Color', [0.9290 0.6940 0.1250])
xlabel('t', 'fontsize', 12)
ylabel('y(t)', 'fontsize', 12)

nexttile
plot(t, u(3, :), 'Color', [0.4660 0.6740 0.1880])
xlabel('t', 'fontsize', 12)
ylabel('z(t)', 'fontsize', 12)

xq = (100 * round(min(maxima)) : 100 * round(max(maxima))) / 100;

figure(2)
hold on
plot(xq, xq, ':', 'Color', 'k') 
scatter(maxima(1 : end - 1), maxima(2 : end),...
    '.', 'MarkerEdgeColor', [0.4940 0.1840 0.5560])
hold off
title(['Lorenz map from t=', num2str(t0),...
    ' to t=', num2str(tf)], 'fontsize', 14)
xlabel('$z_{n}$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$z_{n+1}$', 'interpreter', 'latex', 'fontsize', 14)

figure(3)
scatter3(u(1, :), u(2, :), u(3, :), 0.1, '.')
title(['Lorenz equations attractor from t=', num2str(t0),...
    ' to t=', num2str(tf)], 'fontsize', 14)
xlabel('x(t)', 'fontsize', 12)
ylabel('y(t)', 'fontsize', 12)
zlabel('z(t)', 'fontsize', 12)

cobwebdomain = zeros(1, 2 * (length(maxima) - 1));
cobwebrange = zeros(1, 2 * (length(maxima) - 1));
for i = 1 : length(maxima) - 1
    cobwebdomain(2 * i - 1 : 2 * i) = maxima(i : i + 1);
    cobwebrange(2 * i - 1 : 2 * i) = maxima(i + 1);
end

figure(4)
hold on 
plot(xq, xq, ':', 'Color', 'k') 
plot(cobwebdomain, cobwebrange, 'LineWidth', 0.1) 
scatter(maxima(1 : end - 1), maxima(2 : end),...
    '.', 'MarkerEdgeColor', [0.4940 0.1840 0.5560])
hold off
title(['Lorenz map from t=', num2str(t0),...
    ' to t=', num2str(tf), ' (w/ cobwebs)'], 'fontsize', 14)
xlabel('$z_{n}$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$z_{n+1}$', 'interpreter', 'latex', 'fontsize', 14)