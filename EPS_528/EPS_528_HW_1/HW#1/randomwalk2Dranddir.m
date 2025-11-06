function [ ] = randomwalk2Dranddir(N, I, T, t, X, Y, D, R, reflecting)
close all

cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\randomwalks\2D'
%N = number of agents
%I = number of infected at the start
%T = number of timesteps
%t = recovery time 
%X = horizontal length of domain
%Y = vertical length of domain
%D = expected value of distance moving at each step
%R = infection radius
%BCs are periodic by default. However, if `reflecting` is set to `true` or
%`1`, then we will have reflecting BCs.

%start with agents in random locations
loc = [X; Y] .* rand(2, N);
infections = zeros(1, N); %to track infections
%zero is healthy, nonzero is infected
infections(randperm(N, I)) = t; %set I random agent(s) to be infected
numberinfected = zeros(1, T); %to track number infected
numberinfected(1) = sum(infections > 0);

moviename = ['randomwalk2Dranddir_', num2str(N), '_', num2str(I), '_',...
    num2str(X), 'x', num2str(Y), '_', num2str(D), '_', num2str(R), '_', ...
    num2str(reflecting), '.mp4'];

size = 10;

figure(1)
vidfile = VideoWriter(moviename, 'MPEG-4');
vidfile.FrameRate = 3;
open(vidfile);
scatter(loc(1, :), loc(2, :), size, ceil(infections / t), 'filled')
colormap cool
axis equal
title(['Day ', num2str(1), ' (Percent infected: ', num2str(100 * numberinfected(1) / N), '%)'])
xlim([0 X])    
ylim([0 Y])
xlabel('x')
ylabel('y')
colorbar('Ticks', [0, 1], 'TickLabels', {'Susceptible', 'Infected'})
pause
for i = 2 : T
    infections = max(0, infections - 1); %allow infection to pass
    theta = 2 * pi * rand(1, N);
    Drand = D * sqrt(pi) * erfinv(rand(1, N));
    if reflecting
        loc(1, :) = abs(mod(loc(1, :) + Drand .* cos(theta) - X, 2 * X) - X);
        loc(2, :) = abs(mod(loc(2, :) + Drand .* sin(theta) - Y, 2 * Y) - Y);
    else
        loc(1, :) = mod(loc(1, :) + Drand .* cos(theta), X);
        loc(2, :) = mod(loc(2, :) + Drand .* sin(theta), Y);
    end
    for j = find(infections > 0) % go through infected agents
        distances = sqrt((loc(1, :) - loc(1, j)) .^ 2 + (loc(2, :) -...
            loc(2, j)) .^ 2);
        % check which were infected (but remove already infected)
        infections((distances < R) & (infections == 0)) = t; 
    end
    numberinfected(i) = sum(infections > 0);
    
    scatter(loc(1, :), loc(2, :), size, ceil(infections / t), 'filled')
    axis equal
    title(['Day ', num2str(i), ' (Percent infected: ', num2str(100 * numberinfected(i) / N), '%)'])
    xlim([0 X])    
    ylim([0 Y])
    xlabel('x')
    ylabel('y')
    colorbar('Ticks', [0, 1], 'TickLabels', {'Susceptible', 'Infected'})
    
    pause(0.5)
    drawnow
    F = getframe(gcf); 
    writeVideo(vidfile, F);
end
close(vidfile)

infectedmean = mean(numberinfected);

figure(2)
plot(1 : T, numberinfected, 0 : T, infectedmean * ones(1, T + 1))
title(['Infected vs. time (Mean infected: ', num2str(infectedmean), ')'])
xlabel('Time [days]')
ylabel('Number infected')

cd ..
