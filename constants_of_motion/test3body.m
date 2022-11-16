tic
M = [1 0.9 1];
G = 1 ;
u0 = [1 -0.3 1 1 2 0 0 0 0 0 0 0];
ti = 0;
tf = 100;
tol = 1e-10;
sketch = 0;
snrf = 250;

labels = {'Hamiltonian', 'linear momentum (x)',...
    'linear momentum (y)', 'angular momentum'};

snrq = 10 : 10 : snrf;
N = length(snrq);
sigmaerrvq = zeros(4, N);
choiceerrvq = zeros(4, N);
sigmaerrdq = zeros(4, N);
choiceerrdq = zeros(4, N);
for j = 1 : 4
    % get perfect data
    [~, choicev] = variance3body(1, u0,...
    ti, tf, M, G, tol, Inf, 0, j);
    [~, choiced] = deviation3body(1, u0,...
    ti, tf, M, G, tol, Inf, 0, j);
    for i = 1 : length(snrq)
        [sigmavi, choicevi] = variance3body(1, u0, ti, tf,...
            M, G, tol, snrq(i), sketch, j);
        sigmaerrvq(j, i) = abs(sigmavi);
        choiceerrvq(j, i) = min(norm(choicev - choicevi), ...
            norm(choicev + choicevi));
        [sigmadi, choicedi] = deviation3body(1, u0, ti, tf,...
            M, G, tol, snrq(i), sketch, j);
        sigmaerrdq(j, i) = abs(sigmadi);
        choiceerrdq(j, i) = min(norm(choiced - choicedi), ...
            norm(choiced + choicedi));
        disp(['Finished case for signal-to-noise ratio ',...
            num2str(snrq(i)), '.'])
    end
    disp(['Finished ', num2str(labels{j}), ' cases.'])
end

for i = 1 : 4
    figure(2 * i - 1)
    plot(snrq, [sigmaerrvq(i, :); sigmaerrdq(i, :)])
    legend('Variance functional', 'Deviation functional')
    title(['Weighted norm for ', labels{i}])

    figure(2 * i)
    plot(snrq, [choiceerrvq(i, :); choiceerrdq(i, :)])
    legend('Variance functional', 'Deviation functional')
    title(['Error in coefficients for ', labels{i}])
end
toc