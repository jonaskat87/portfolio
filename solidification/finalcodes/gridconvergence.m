dxq = 0.008 : -0.001 : 0.004;
N = length(dxq);
phiL2errors = zeros(1, N);
TL2errors = zeros(1, N); 
timeerrors = zeros(1, N); 
heighterrors = zeros(1, N);

d = 1;
scaling = 4;
Tend = 800 + 273.15;
dtq = 60 * 60;

load('BIGemplacementresults.mat')
phiuber = all{2};
Tuber = all{3};
datalastuber = all{11};
lasttime = mean(datalastuber(1));
lastheight = mean(datalastuber(2));
dxuber = 4e-4;
phiweight = sqrt(dxuber) * norm(phiuber, 'fro');
Tweight = sqrt(dxuber) * norm(Tuber, 'fro');

for i = 1 : N
    all0 = simulation2fortestingtemp(scaling, dxq(i), d, 2, 0, dtq,...
 dtq, Tend, false);
    phi0 = all0{2};
    T0 = all0{3};
    datalast0 = all0{11};
    lasttime0 = mean(datalast0(1));
    lastheight0 = mean(datalast0(2));
    factor = dxq(i) / dxuber;
    if size(phi0) == size(phiuber(factor : factor : end, 1 : length(phi0(1, :))))
        phiL2errors(i) = sqrt(dxq(i)) * norm(phi0 - phiuber(factor...
            : factor : end, 1 : length(phi0(1, :))), 'fro') / phiweight;
        TL2errors(i) = sqrt(dxq(i)) * norm(T0 - Tuber(factor : factor...
            : end, 1 : length(phi0(1, :))), 'fro') / Tweight;
    else
        for j = 1 : factor
            if j == factor
               disp('Something wrong!')
               return 
            end
            if size(phi0) == size(phiuber(factor - j : factor : end, 1 : length(phi0(1, :))))
                phiL2errors(i) = sqrt(dxq(i)) * norm(phi0 - phiuber(factor...
                    - j : factor : end, 1 : length(phi0(1, :))), 'fro') / phiweight;
                TL2errors(i) = sqrt(dxq(i)) * norm(T0 - Tuber(factor -...
                    j : factor : end, 1 : length(phi0(1, :))), 'fro') / Tweight;   
                break
            elseif size(phi0) == size(phiuber(factor : factor : end - j, 1 : length(phi0(1, :))))
                phiL2errors(i) = sqrt(dxq(i)) * norm(phi0 - phiuber(factor...
                    : factor : end - j, 1 : length(phi0(1, :))), 'fro') / phiweight;
                TL2errors(i) = sqrt(dxq(i)) * norm(T0 - Tuber(factor...
                    : factor : end - j, 1 : length(phi0(1, :))), 'fro') / Tweight;   
            end
        end
    end
    timeerrors(i) = abs(lasttime - lasttime0) / abs(lasttime);
    heighterrors(i) = abs(lastheight - lastheight0) / abs(lastheight);
    disp(['I just finished the case for dx=', num2str(dxq(i)), '.']) 
end

figure
tiledlayout(4,1)
nexttile
plot(dxq, phiL2errors)
grid on
title('Phi errors')
xlabel('$\Delta x$', 'interpreter', 'latex')
ylabel('Relative L2 error')

nexttile
plot(dxq, TL2errors)
grid on
title('T errors')
xlabel('$\Delta x$', 'interpreter', 'latex')
ylabel('Relative L2 error')

nexttile
plot(dxq, timeerrors)
grid on
title('Solidification time errors')
xlabel('$\Delta x$', 'interpreter', 'latex')
ylabel('Relative error')

nexttile
plot(dxq, heighterrors)
grid on
title('Solidification height errors')
xlabel('$\Delta x$', 'interpreter', 'latex')
ylabel('Relative error')

pphi = log(phiL2errors(2 : end) ./ phiL2errors(1 : end - 1)) ./ (log(dxq(2 : end) ./ dxq(1 : end - 1)));
disp(['Rate for phi errors: ', num2str(pphi), '.'])

pT = log(TL2errors(2 : end) ./ TL2errors(1 : end - 1)) ./ (log(dxq(2 : end) ./ dxq(1 : end - 1)));
disp(['Rate for T errors: ', num2str(pT), '.'])

ptime = log(timeerrors(2 : end) ./ timeerrors(1 : end - 1)) ./ (log(dxq(2 : end) ./ dxq(1 : end - 1)));
disp(['Rate for solidification time errors: ', num2str(ptime), '.'])

pheight = log(heighterrors(2 : end) ./ heighterrors(1 : end - 1)) ./ (log(dxq(2 : end) ./ dxq(1 : end - 1)));
disp(['Rate for solidification height errors: ', num2str(pheight), '.'])