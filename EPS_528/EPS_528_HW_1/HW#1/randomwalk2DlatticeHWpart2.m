function [ ] = randomwalk2DlatticeHWpart2(N, M, Q, L, n)
numsimulations = 10;
data = zeros(numsimulations, n);

for i = 1 : numsimulations
    data(i, :) = randomwalk2DlatticeHWpart1(N, M, Q, L, n);
    disp(i)
end

numberinfected = mean(data); 

figure
plot(1 : n, numberinfected)
title(['Infected vs. time, M=', num2str(M)])
xlabel('Time [days]')
ylabel('Number infected')