function [ ] = simulation2scatter(N, Tend, actual)

% mq = [0.1 : 0.1 : 1 2 3 4 5 10 15 20];
% range = 1 : length(mq);
% tsolidfull = zeros(2, length(mq) * N);
% tsolid1st = zeros(2, length(mq) * N);
% zfull = zeros(2, length(mq) * N);
% z1st = zeros(2, length(mq) * N);
% for i = range
%     datai = simulation2extraction(N, i, Tend, actual);
%     tsolidfull(:, N * (i - 1) + 1 : N * i) = [datai(1, :);...
%         datai(2, :) - datai(1, :)];
%     tsolid1st(:, N * (i - 1) + 1 : N * i) = [datai(1, :); datai(3, :)];
%     zfull(:, N * (i - 1) + 1 : N * i) = [datai(1, :); datai(5, :)];
%     z1st(:, N * (i - 1) + 1 : N * i) = [datai(1, :); datai(6, :)];
%     disp(['Finished h=', num2str(mq(i)), 'm.'])
% end
% 
% save scatterdata.mat tsolidfull tsolid1st zfull z1st
% return

load scatterdata.mat tsolidfull tsolid1st zfull z1st

earlierrange = N + 1 : 11 * N;
laterrange = 11 * N + 1 : length(tsolidfull(1, :));
% range = 1 : length(tsolidfull(1, :));

tol = 1e-15;
F1 = @(x, xdata) x(1) - x(2) * (1 - exp(-x(3) * xdata)) -...
    x(4) * (1 - exp(-x(5) * xdata .^ 2));
x0 = [4 3 20 2 20];
lb = [0 0 0 0 0];
options = optimoptions('lsqcurvefit', 'FunctionTolerance', tol,...
    'OptimalityTolerance', tol, 'MaxIterations', 5e3,...
    'MaxFunctionEvaluations', 5e3, 'StepTolerance', tol);
[xf1earlier, resnorm1earlier, ~, ~, ~] = lsqcurvefit(F1, x0, tsolidfull(1, earlierrange), tsolidfull(2, earlierrange), lb, [], options);
disp(xf1earlier)
[xf1later, resnorm1later, ~, ~, ~] = lsqcurvefit(F1, x0, tsolidfull(1, laterrange), tsolidfull(2, laterrange), lb, [], options);
disp(xf1later)

L2error1earlier = sqrt(resnorm1earlier) / sqrt(sum(tsolidfull(2, earlierrange) .^ 2));
disp(L2error1earlier)
L2error1later = sqrt(resnorm1later) / sqrt(sum(tsolidfull(2, laterrange) .^ 2));
disp(L2error1later)

F2 = @(x, xdata) x(1) - x(2) * (1 - exp(-x(3) * xdata)) -...
    (x(1) - x(2) - 1) * (1 - exp(-x(4) * xdata .^ 2));
x0 = [4 3 20 20];
lb = [0 0 0 0];
[xf2earlier, resnorm2earlier, ~, ~, ~] = lsqcurvefit(F2, x0, tsolid1st(1, earlierrange), tsolid1st(2, earlierrange), lb, [], options);
disp(xf2earlier)
[xf2later, resnorm2later, ~, ~, ~] = lsqcurvefit(F2, x0, tsolid1st(1, laterrange), tsolid1st(2, laterrange), lb, [], options);
disp(xf2later)

L2error2earlier = sqrt(resnorm2earlier) / sqrt(sum(tsolid1st(2, earlierrange) .^ 2));
disp(L2error2earlier)
L2error2later = sqrt(resnorm2later) / sqrt(sum(tsolid1st(2, laterrange) .^ 2));
disp(L2error2later)

F3 = @(x, xdata) x(1) + x(2) * tanh(x(3) + x(4) * log(xdata));
x0 = [4 3 2 1];
lb = [0 0 0 0];
options = optimoptions('lsqcurvefit', 'FunctionTolerance', tol,...
    'OptimalityTolerance', tol, 'MaxIterations', 5e3,...
    'MaxFunctionEvaluations', 5e3, 'StepTolerance', tol);
[xf3earlier, resnorm3earlier, ~, ~, ~] = lsqcurvefit(F3, x0, zfull(1, earlierrange), zfull(2, earlierrange), lb, [], options);
disp(xf3earlier)
[xf3later, resnorm3later, ~, ~, ~] = lsqcurvefit(F3, x0, zfull(1, laterrange), zfull(2, laterrange), lb, [], options);
disp(xf3later)

L2error3earlier = sqrt(resnorm3earlier) / sqrt(sum(zfull(2, earlierrange) .^ 2));
disp(L2error3earlier)
L2error3later = sqrt(resnorm3later) / sqrt(sum(zfull(2, laterrange) .^ 2));
disp(L2error3later)

earlierlaterrange = zeros(1, 7 * length(earlierrange) / N);
laterlaterrange = zeros(1, 8 * length(laterrange) / N);
for i = 2 : N
    if i == 2
    	laterlaterrange(6 * (i - 2) + (1 : 6)) = laterrange(i : N : length(laterrange));
    else
        earlierlaterrange(10 * (i - 3) + (1 : 10)) = earlierrange(i : N : length(earlierrange));
        laterlaterrange(6 * (i - 2) + (1 : 6)) = laterrange(i : N : length(laterrange));
    end
end

F4 = @(x, xdata) x(1) - x(2) * tanh(x(3) + x(4) * log(xdata));
x0 = [4 3 2 1];
lb = [0 0 0 0];
options = optimoptions('lsqcurvefit', 'FunctionTolerance', tol,...
    'OptimalityTolerance', tol, 'MaxIterations', 5e3,...
    'MaxFunctionEvaluations', 5e3, 'StepTolerance', tol);
[xf4earlier, resnorm4earlier, ~, ~, ~] = lsqcurvefit(F4, x0, z1st(1, earlierlaterrange), z1st(2, earlierlaterrange), lb, [], options);
disp(xf4earlier)
[xf4later, resnorm4later, ~, ~, ~] = lsqcurvefit(F4, x0, z1st(1, laterlaterrange), z1st(2, laterlaterrange), lb, [], options);
disp(xf4later)

L2error4earlier = sqrt(resnorm4earlier) / sqrt(sum(z1st(2, earlierlaterrange) .^ 2));
disp(L2error4earlier)
L2error4later = sqrt(resnorm4later) / sqrt(sum(z1st(2, laterlaterrange) .^ 2));
disp(L2error4later)

figure(1)
tiledlayout(2, 1)
nexttile
grid on
hold on
scatter(tsolidfull(1, earlierrange), tsolidfull(2, earlierrange), 'r')
scatter(tsolidfull(1, laterrange), tsolidfull(2, laterrange), 'b')
plot(0 : 0.01 : 18, F1(xf1earlier, 0 : 0.01 : 18), 'r')
plot(0 : 0.01 : 18, F1(xf1later, 0 : 0.01 : 18), 'b')
hold off
set(gca, 'Xscale', 'log')
title('Time to solidification across domain')
xlabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
ylabel('$(t_{solidification}-t_{emp})/t_{h}$', 'interpreter', 'latex')
legend('Smaller', 'Larger', ['Fit for smaller sizes (relative L2 error: ',...
    num2str(L2error1earlier), ')'],...
    ['Fit for larger sizes (relative L2 error: ', num2str(L2error1later), ')'])

nexttile
grid on
hold on
scatter(tsolid1st(1, earlierrange), tsolid1st(2, earlierrange), 'r')
scatter(tsolid1st(1, laterrange), tsolid1st(2, laterrange), 'b')
plot(0 : 0.01 : 18, F2(xf2earlier, 0 : 0.01 : 18), 'r')
plot(0 : 0.01 : 18, F2(xf2later, 0 : 0.01 : 18), 'b')
hold off
set(gca, 'Xscale', 'log')
title('Time to solidification for first lobe only')
xlabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
ylabel('$t_{solidification}/t_{h}$', 'interpreter', 'latex')
legend('Smaller', 'Larger', ['Fit for smaller sizes (relative L2 error: ',...
    num2str(L2error2earlier), ')'],...
    ['Fit for larger sizes (relative L2 error: ', num2str(L2error2later), ')'])

figure(2)
tiledlayout(2, 1)
nexttile
grid on
hold on
scatter(zfull(1, earlierrange), zfull(2, earlierrange), 'r')
scatter(zfull(1, laterrange), zfull(2, laterrange), 'b')
plot(0 : 0.01 : 18, F3(xf3earlier, 0 : 0.01 : 18), 'r')
plot(0 : 0.01 : 18, F3(xf3later, 0 : 0.01 : 18), 'b')
hold off
set(gca, 'Xscale', 'log')
title('Height where solidification completed across domain')
xlabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
ylabel('$z/h$', 'interpreter', 'latex')
legend('Smaller', 'Larger', ['Fit for smaller sizes (relative L2 error: ',...
    num2str(L2error3earlier), ')'],...
    ['Fit for larger sizes (relative L2 error: ', num2str(L2error3later), ')'],...
    'Location', 'southeast')

nexttile
grid on
hold on
h1 = scatter(z1st(1, earlierrange), z1st(2, earlierrange), 'r');
h2 = scatter(z1st(1, laterrange), z1st(2, laterrange), 'b');
h3 = plot(0 : 0.01 : 0.13, F3(xf3earlier, 0 : 0.01 : 0.13), 'r');
plot([0.13 0.26 0.53], [F3(xf3earlier, 0.13)...
    F4(xf4earlier, 0.26) F4(xf4earlier, 0.53)], '--', 'color', 'r')
plot(0.53 : 0.01 : 18, F4(xf4earlier, 0.53 : 0.01 : 18), 'r')
h4 = plot(0 : 0.01 : 0.13, F3(xf3later, 0 : 0.01 : 0.13), 'b');
plot(0.26 : 0.01 : 18, F4(xf4later, 0.26 : 0.01 : 18), 'b')
plot([0.13 0.26], [F3(xf3later, 0.13) F4(xf4later, 0.26)], '--', 'color', 'b')
hold off
set(gca, 'Xscale', 'log')
title('Height where solidification completed in the first lobe')
xlabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
ylabel('$z/h$', 'interpreter', 'latex')
legend([h1 h2 h3 h4], {'Smaller', 'Larger',...
    ['Fit for smaller sizes (relative L2 error: ', num2str(L2error4earlier),...
    ')'], ['Fit for larger sizes (relative L2 error: ',...
    num2str(L2error4later), ')']}, 'Location', 'east')