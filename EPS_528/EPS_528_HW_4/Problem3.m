cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\EPS_528_HW_4'

N = 5;

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

xq = (100 * round(min(maxima)) : 100 * round(max(maxima))) / 100;

figure(1)
tiledlayout(N, 1)
for i = 1 : N
    nexttile
    hold on
    plot(xq, xq, ':', 'Color', 'k') 
    scatter(maxima(1 : end - i), maxima(i + 1 : end),...
        '.', 'MarkerEdgeColor', [0.4940 0.1840 0.5560])
    hold off
    if i == 1
        title(['Lorenz map from t=', num2str(t0),...
            ' to t=', num2str(tf)], 'fontsize', 14)
    end
    xlabel('$z_{n}$', 'interpreter', 'latex', 'fontsize', 14)
    ylabel(['$z_{n+', num2str(i), '}$'],...
        'interpreter', 'latex', 'fontsize', 14)
end
