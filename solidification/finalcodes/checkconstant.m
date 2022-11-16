% cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\newemplacementresults'
% d = 1;
% scaling = 5;
% dx = round(sqrt(10) ^ (log10(0.1 / 2) - 5), 3);
% Tend = 1000 + 273.15;
% alpha = 3.75e-7;
% dtq = 60;
% % dtq = 60 * 10;
% 
% case1 = simulation2fortestingtemp(scaling, dx, d, 0.1, 0, dtq, dtq, Tend, false);
% fname = 'propconstant.mat';
% % fname = 'propconstantmany.mat';
% save(fname, 'case1');
% m = matfile(fname, 'Writable', true);
% disp(['I just finished the case for a ', num2str(0.1), ' meter lobe.'])
% 
% mq = [0.1 : 0.1 : 1 2 3 4 5 10 15 20];
% % mq = [0.1 : 0.05 : 0.5 0.6 0.8 1 2 3 4 5 7 10];
% for i = 2 : length(mq)
%     dx = min(round(sqrt(10) ^ (log10(mq(i)) - 5), 3), 0.01);
%     dtq = 60 * (mq(i) / 0.1) ^ 2;
%     all = simulation2fortestingtemp(scaling, dx, d, mq(i), 0, dtq, dtq, Tend, false);
%     m.(['case', num2str(i)]) = all;
%     disp(['I just finished the case for a ', num2str(mq(i)), ' meter lobe.'])  
% end
% return

% cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\newemplacementresults'
% 
% % mq = [0.1 : 0.1 : 1 2 3 4 5 10 15 20];
% mq = [0.1 : 0.05 : 0.5 0.6 0.8 1 2 3 4 5 7 10];
% 
% alpha = 3.75e-7;
% 
% % fname = 'propconstant.mat';
% fname = 'propconstantmany.mat';
% m = matfile(fname);
% 
% n = length(mq);
% data = zeros(1, n);
% for i = 1 : length(mq)
%     all = m.(['case', num2str(i)]);
%     datalastoverall = all{11};
%     th = mq(i) ^ 2  / alpha;
%     data(i) = mean(datalastoverall(1, :)) / th;
% end
% 
% save propconstantmanydata.mat data
% return

load propconstantmanydata.mat data

mq = [0.1 : 0.05 : 0.5 0.6 0.8 1 2 3 4 5 7 10];
datamean = mean(data);
datastd = std(data);

hinterp = 0.1 : 0.05 : 10;
datainterp = interp1(mq, data, hinterp, 'pchip');
% figure
% plot(hinterp, 0.006 * hinterp .^ (-1 / 2) + 0.037, hinterp, datainterp)
% legend

tol = 1e-15;
% F = @(x, xdata) x(1) * (xdata .^ (x(2))) + x(3);
% x0 = [0.006 -0.5 0.037];
% lb = [0 -Inf 0];
% ub = [Inf 0 Inf];
F = @(x, xdata) x(1) * (xdata .^ (x(2))) + x(3) + x(4) * exp(x(5) * xdata);
x0 = [0.006 -0.5 0.037 0.05 -10];
lb = [0 -Inf 0 0 -Inf];
ub = [Inf 0 Inf Inf 0];
options = optimoptions('lsqcurvefit', 'FunctionTolerance', tol, 'OptimalityTolerance', tol, 'MaxFunctionEvaluations', 5e3);
[xf, resnorm, res, exitflag, output] = lsqcurvefit(F, x0, hinterp, datainterp, lb, ub, options);
disp(xf)

L2error = sqrt(sum((datainterp - F(xf, hinterp)) .^ 2)) /...
    sqrt(sum(datainterp .^ 2));

Fapprox = @(x, xdata) x(1) * (xdata .^ (x(2))) + x(3);
% 1 percent relative error: 0.26344

figure
loglog(0.1 : 0.01 : 20, (F(xf, 0.1 : 0.01 : 20) - Fapprox(xf, 0.1 : 0.01 : 20)) ./ F(xf, 0.1 : 0.01 : 20))

figure
% plot(mq, data, 'r*', mq, datamean * ones(1, n), 'b')
patch([0.27 20 20 0.27], [F(xf, 20) F(xf, 20) F(xf, 0.1) + 0.005 F(xf, 0.1) + 0.005], [255 255 102] / 255, 'EdgeColor', [255 255 102] / 255)
hold on
patch([0.1 0.27 0.27 0.1], [F(xf, 20) F(xf, 20) F(xf, 0.1) + 0.005 F(xf, 0.1) + 0.005], [255 153 51] / 255, 'EdgeColor', [255 153 51] / 255)
loglog(mq, data, 'r*', 0.1 : 0.01 : 20, F(xf, 0.1 : 0.01 : 20), 'b-','Linewidth',1.5)
x1 = xline(0.27, '--k', 'Weakly nonlinear');
x1.LabelVerticalAlignment = 'top';
x2 = xline(0.27, '--k', 'Strongly nonlinear', 'LabelHorizontalAlignment', 'left');
x2.LabelVerticalAlignment = 'top';
hold off
set(gca, 'Xscale', 'log', 'Yscale', 'log')
title('Proportionality constant', 'fontsize', 16)
xlabel('$h$ (meters)', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$t_{h}\alpha/h^{2}$', 'interpreter', 'latex', 'fontsize', 14)
legend('Exact', ['Fit (relative L2 error: ', num2str(round(L2error,3)), ')'])
% text(2, 0.0485, ['standard deviation: ', num2str(datastd)])

figure
% plot(mq, data, 'r*', mq, datamean * ones(1, n), 'b')
loglog(0.1 : 0.1 : 100, F(xf, 0.1 : 0.1 : 100), 'b-')
title('Fit extrapolation', 'fontsize', 16)
xlabel('$h$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$t_{h}\alpha/h^{2}$', 'interpreter', 'latex', 'fontsize', 14)