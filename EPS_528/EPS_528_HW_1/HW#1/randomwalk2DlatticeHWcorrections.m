function [ ] = randomwalk2DlatticeHWcorrections(N, M, Q, L, n)
numsimulations = 10;
data = zeros(numsimulations, n);

for i = 1 : numsimulations
    data(i, :) = randomwalk2DlatticeHWpart1(N, M, Q, L, n);
    disp(i)
end

func1 = @(x) x * (1 + (4 / 5) * ((M - x) / (N ^ 2)));
func2 = @(x) ((L - 1) / L) * x * (1 + (4 / 5) *...
    ((M - x) / (N ^ 2)) + x / (L * N ^ 2));
I0 = Q;
I1 = Q;
Idifference = zeros(2, n);
Idifference(1) = Q;
for i = 2 : n
    if i < L
        I1 = func1(I1);
    else
    	I1 = func2(I1);
    end
    I0 = func2(I0);
    Idifference(1, i) = I0;
    Idifference(2, i) = I1;
end

numberinfected = mean(data); 

figure
plot(1 : n, numberinfected, 1 : n, Idifference(1, :),...
    1 : n, Idifference(2, :))
title(['Infected vs. time, M=', num2str(M)])
xlabel('Time [days]')
ylabel('Number infected')
legend('Agent-based model', 'Discrete model 0', 'Discrete model 1')