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
[maxdomain, I] = sort(maxima(1 : end - 1));
maxvalue = maxima(2 : end);
maxvalue = maxvalue(I);

fprimefirst = (maxvalue(3 : end) - maxvalue(2 : end - 1))...
    ./ (maxdomain(3 : end) - maxdomain(2 : end - 1));
fprimesecond = (maxvalue(2 : end - 1) - maxvalue(1 : end - 2))...
    ./ (maxdomain(2 : end - 1) - maxdomain(1 : end - 2));
fprime = 0.5 * (fprimefirst + fprimesecond);
x = maxdomain(2 : end - 1);
y = abs(fprime);
threshold = 10;

xq = (100 * round(min(x)) :...
    100 * round(max(x))) / 100;

figure(1)
hold on
scatter(x(y < threshold), y(y < threshold), '.',...
    'MarkerEdgeColor', [0.4940 0.1840 0.5560])
yline(1, '--k', '$\left|f^{\prime}(z)\right|=1$', 'interpreter', 'latex');
hold off
title(['$\left|f^{\prime}(z)\right|$', ' from t=',...
    num2str(t0), ' to t=', num2str(tf), ' with ',...
    '$\left|f^{\prime}(z)\right|<$', ' ', num2str(threshold)],...
    'interpreter', 'latex', 'fontsize', 14)
xlabel('$z$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$\left|f^{\prime}(z)\right|$',...
    'interpreter', 'latex', 'fontsize', 14)