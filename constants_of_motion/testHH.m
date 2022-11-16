tic
u0 = [0.1, 0.4, -0.3, -0.2];
ti = 0;
tf = 100;
lambda = 0.9;
tol = 1e-12;
sketch = 0;
snrf = 250;

% get perfect data
[sigmav, choicev] = varianceHH(1, u0, ti, tf, lambda, tol, Inf, 0);
[sigmad, choiced] = deviationHH(1, u0, ti, tf, lambda, tol, Inf, 0);

snrq = 10 : 1 : snrf;
N = length(snrq);
sigmaerrvq = zeros(1, N);
choiceerrvq = zeros(1, N);
sigmaerrdq = zeros(1, N);
choiceerrdq = zeros(1, N);
for i = 1 : length(snrq)
    [sigmavi, choicevi] = varianceHH(1, u0, ti, tf,...
        lambda, tol, snrq(i), sketch);
    sigmaerrvq(i) = abs(sigmavi - sigmav);
    choiceerrvq(i) = min(norm(choicev - choicevi),...
        norm(choicev + choicevi));
    [sigmadi, choicedi] = deviationHH(1, u0, ti, tf,...
        lambda, tol, snrq(i), sketch);
    sigmaerrdq(i) = abs(sigmadi - sigmad);
    choiceerrdq(i) = min(norm(choiced - choicedi), ...
        norm(choiced + choicedi));
    disp(['Finished case for signal-to-noise ratio ',...
        num2str(snrq(i)), '.'])
end

figure(3)
plot(snrq, [sigmaerrvq; sigmaerrdq])
legend('Variance functional', 'Deviation functional')
title('Error in norm')

figure(4)
plot(snrq, [choiceerrvq; choiceerrdq])
legend('Variance functional', 'Deviation functional')
title('Error in coefficients')
toc