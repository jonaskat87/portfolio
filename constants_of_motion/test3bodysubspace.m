tic
M = [1 0.9 1];
G = 1 ;
u0 = [1 -0.3 1 1 2 0 0 0 0 0 0 0];
ti = 0;
tf = 100;
tol = 1e-12;
sketch = 0;
snrf = 250;

labels = {'Hamiltonian', 'linear momentum (x)',...
    'linear momentum (y)', 'angular momentum'};

snrq = 10 : 2 : snrf;
N = length(snrq);
sigmaerrvq = zeros(1, N);
choiceerrvq = zeros(1, N);
sigmaerrdq = zeros(1, N);
choiceerrdq = zeros(1, N);
% get perfect data
[~, choicev] = variance3body(4, u0,...
    ti, tf, M, G, tol, Inf, 0, 5);
Pv = choicev * ((transpose(choicev) * choicev) \ transpose(choicev));
[~, choiced] = deviation3body(4, u0,...
    ti, tf, M, G, tol, Inf, 0, 5);
Pd = choiced * ((transpose(choiced) * choiced) \ transpose(choiced));
for i = 1 : length(snrq)
    [sigmavi, choicevi] = variance3body(4, u0, ti, tf,...
        M, G, tol, snrq(i), sketch, 5);
    sigmaerrvq(i) = norm(sigmavi);
    Pvi = choicevi * ((transpose(choicevi) *...
        choicevi) \ transpose(choicevi));
    choiceerrvq(i) = norm(Pv - Pvi);
    [sigmadi, choicedi] = deviation3body(4, u0, ti, tf,...
        M, G, tol, snrq(i), sketch, 5);
    sigmaerrdq(i) = norm(sigmadi);
    Pdi = choicedi * ((transpose(choicedi) *...
        choicedi) \ transpose(choicedi));    
    choiceerrdq(i) = norm(Pd - Pdi);
    disp(['Finished case for signal-to-noise ratio ',...
        num2str(snrq(i)), '.'])
end

figure(1)
plot(snrq, [sigmaerrvq; sigmaerrdq])
legend('Variance functional', 'Deviation functional')
title('Weighted norm')

figure(2)
plot(snrq, [choiceerrvq; choiceerrdq])
legend('Variance functional', 'Deviation functional')
title('Error in coefficients (Grassmannian distance)')
toc