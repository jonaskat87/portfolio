function [ ] = randomwalk2Dlattice(N, M, Q, L, n)
close all

cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\randomwalks\2D'
%N = linear dimension (lattice is size N^2)
%M = total population number
%Q = number of initial infected people
%L = days to recover
%n = number of total iterations

%start with agents in random locations
loc = randi(N, 2, M);
infections = zeros(1, M); %to track infections
%zero is healthy, nonzero is infected
infections(randperm(M, Q)) = L; %set Q random agent(s) to be infected
numberinfected = zeros(1, n); %to track number infected
numberinfected(1) = sum(infections > 0);

%sample directions (to add on at each iteration)
dirq = [0 0 0 1 -1; 0 1 -1 0 0];

moviename = ['randomwalk2Dlattice_', num2str(N), '_', num2str(M), '_',...
    num2str(Q), '_', num2str(L), '_', num2str(n), '.mp4'];

size = 10;

figure(1)
vidfile = VideoWriter(moviename, 'MPEG-4');
vidfile.FrameRate = 3;
open(vidfile);
scatter(loc(1, :), loc(2, :), size, ceil(infections / L), 'filled')
colormap cool
axis equal
title(['Day ', num2str(1), ' (Percent infected: ', num2str(100 * numberinfected(1) / M), '%)'])
xlim([0 N])    
ylim([0 N])
xlabel('x')
ylabel('y')
colorbar('Ticks', [0, 1], 'TickLabels', {'Susceptible', 'Infected'})
grid on
pause
for i = 2 : n
    infections = max(0, infections - 1); %allow infection to pass
    loc = mod(loc + dirq(:, randi(5, 1, M)), N);
    for j = find(infections > 0) % go through infected agents
        infloc = loc(:, j);
        % check where locations line up with infected
        equality = (loc == infloc); 
        % check which were infected (but remove already infected)
        infections((sum(equality) == 2) & (infections == 0)) = L; 
    end
    numberinfected(i) = sum(infections > 0);
    
    scatter(loc(1, :), loc(2, :), size, ceil(infections / L), 'filled')
    axis equal
    title(['Day ', num2str(i), ' (Percent infected: ', num2str(100 * numberinfected(i) / M), '%)'])
    xlim([0 N])    
    ylim([0 N])
    xlabel('x')
    ylabel('y')
    colorbar('Ticks', [0, 1], 'TickLabels', {'Susceptible', 'Infected'})
    grid on
    
    drawnow
    F = getframe(gcf); 
    writeVideo(vidfile, F);
end
close(vidfile)

infectedmean = mean(numberinfected);

figure(2)
plot(1 : n, numberinfected, 0 : n, infectedmean * ones(1, n + 1))
title(['Infected vs. time (Mean infected: ', num2str(infectedmean), ')'])
xlabel('Time [days]')
ylabel('Number infected')

cd ..
