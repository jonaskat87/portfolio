function [ ] = randomwalk2DlatticeHWpart3(N, M, Q, L, n)
numsimulations = 10;
data = zeros(numsimulations, n);

for i = 1 : numsimulations
    data(i, :) = randomwalk2DlatticeHWpart1(N, M, Q, L, n);
    disp(i)
end

func = @(x) ((L - 1) / L) * x * (1 + (4 / 5) *...
    ((M - x) / (N ^ 2)) + x / (L * N ^ 2));
I = Q;
Idifference = zeros(1, n);
Idifference(1) = Q;
for i = 2 : n
    I = func(I);
    Idifference(i) = I;
end

numberinfected = mean(data); 

figure
plot(1 : n, numberinfected, 1 : n, Idifference)
title(['Infected vs. time, M=', num2str(M)])
xlabel('Time [days]')
ylabel('Number infected')
legend('Agent-based model', 'Discrete model 0')