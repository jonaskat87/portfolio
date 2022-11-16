cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\EPS_528_HW_4'

sigma = 10;
b = 8 / 3;
r = 28;
x0 = 0;
y0 = 1;
z0 = 0;
t0 = 0;
tf = 100;
sol0 = lorenz(sigma, r, b, t0, tf, x0, y0, z0);
solp = lorenz(sigma, r, b, t0, tf, x0 + 1e-6, y0, z0);
t0 = sol0.x;
u0 = sol0.y;
tp = solp.x;
up = solp.y;

% create Lorenz map
maxima0 = NaN(1, length(t0));
maximap = NaN(1, length(tp));
samplei0 = u0(3, 1);
sampleip = up(3, 1);
sampleafter0 = u0(3, 2);
sampleafterp = up(3, 2);
for i = 2 : max(length(t0), length(tp)) - 1
    if i < length(t0)
        samplebefore0 = samplei0;
        samplei0 = sampleafter0;
        sampleafter0 = u0(3, i + 1);
        if (samplei0 > samplebefore0) && (samplei0 > sampleafter0)
            maxima0(i) = samplei0;
        end
    end
    if i < length(tp)
        samplebeforep = sampleip;
        sampleip = sampleafterp;    
        sampleafterp = up(3, i + 1);
        if (sampleip > samplebeforep) && (sampleip > sampleafterp)
            maximap(i) = sampleip;
        end
    end
end
maxima0 = maxima0(~(isnan(maxima0)));
maximap = maximap(~(isnan(maximap)));

figure(1)
tiledlayout(3, 1)
nexttile
plot(t0, u0(1, :), tp, up(1, :))
title('Lorenz equation solution', 'fontsize', 14)
xlabel('t', 'fontsize', 12)
ylabel('x(t)', 'fontsize', 12)
legend('original', 'perturbed')

nexttile
plot(t0, u0(2, :), tp, up(2, :))
xlabel('t', 'fontsize', 12)
ylabel('y(t)', 'fontsize', 12)
legend('original', 'perturbed')

nexttile
plot(t0, u0(3, :), tp, up(3, :))
xlabel('t', 'fontsize', 12)
ylabel('z(t)', 'fontsize', 12)
legend('original', 'perturbed')

figure(2)
xlabel('x(t)', 'fontsize', 12)
ylabel('y(t)', 'fontsize', 12)
zlabel('z(t)', 'fontsize', 12)
hold on
scatter3(u0(1, :), u0(2, :), u0(3, :), 0.1, '.')
scatter3(up(1, :), up(2, :), up(3, :), 0.1, '.')
hold off
grid on
title('Lorenz equations attractor', 'fontsize', 14)
legend('original', 'perturbed')

if length(tp) < length(t0)
    u = u0;
    uinterp = zeros(1, length(t0));
    for i = 1 : 3
        uinterp(i, :) = interp1(tp, up(i, :), t0, 'makima');
    end
    tq = t0;
else
    u = up;
    uinterp = zeros(1, length(tp));
    for i = 1 : 3
        uinterp(i, :) = interp1(t0, u0(i, :), tp, 'makima');
    end
    tq = tp;
end
distance = 0;
for i = 1 : 3
    distance = distance + (uinterp(i, :) - u(i, :)) .^ 2;  
end
distance = sqrt(distance);
 
A = [ones(length(tq(round(end / 100) : round(end / 6))),...
    1) tq(round(end / 100) : round(end / 6)).'];
b = log(distance(round(end / 100) : round(end / 6))).';
coeff1 = A \ b;
disp(['The MLE in the first section is ', num2str(coeff1(2))])

A = [ones(length(tq(round(end / 6) : round(35 * end / 100))),...
    1) tq(round(end / 6) : round(35 * end / 100)).'];
b = log(distance(round(end / 6) : round(35 * end / 100))).';
coeff2 = A \ b;
disp(['The MLE in the second section is ', num2str(coeff2(2))])

figure(3)
hold on
plot(tq, log(distance))
plot(tq(round(end / 100) : round(end / 6)), coeff1(1)...
    + coeff1(2) * tq(round(end / 100) : round(end / 6)), '--')
plot(tq(round(end / 6) : round(35 * end / 100)), coeff2(1)...
    + coeff2(2) * tq(round(end / 6) : round(35 * end / 100)), '--')
hold off
title('log(distance) vs. time', 'fontsize', 14)
xlabel('t', 'fontsize', 12, 'interpreter', 'latex')
ylabel(['$\log\left|\left|\mathbf{u}_{\varepsilon}(t)'...
    , '-\mathbf{u}_{0}(t)\right|\right|_{2}$'],...
     'fontsize', 12, 'interpreter', 'latex')
legend('log(distance)', ['slope: ', num2str(coeff1(2))],...
    ['slope: ', num2str(coeff2(2))])
