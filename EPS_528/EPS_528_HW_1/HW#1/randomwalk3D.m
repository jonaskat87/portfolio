function [ ] = randomwalk3D(N, I, T, t, X, Y, Z, D, R, reflecting)
close all

cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\randomwalks\3D'
%N = number of agents
%I = number of infected at the start
%T = number of timesteps
%t = recovery time 
%X, Y, Z are the lengths of the domain in each direction
%D = distance moving at each step
%R = infection radius
%BCs are periodic by default. However, if `reflecting` is set to `true` or
%`1`, then we will have reflecting BCs.

%start with agents in random locations
loc = [X; Y; Z] .* rand(3, N);
infections = zeros(1, N); %to track infections
%zero is healthy, nonzero is infected
infections(randperm(N, I)) = t; %set I random agent(s) to be infected
numberinfected = zeros(1, T); %to track number infected
numberinfected(1) = sum(infections > 0);

moviename = ['randomwalk3D_', num2str(N), '_', num2str(I), '_',...
    num2str(X), 'x', num2str(Y), '_', num2str(D), '_', num2str(R), '_', ...
    num2str(reflecting), '.mp4'];

size = 10;

figure(1)
tile = tiledlayout(1, 2);
title(tile, ['Day ', num2str(1), ' (Percent infected: ', num2str(100 * numberinfected(1) / N), '%)'])
vidfile = VideoWriter(moviename, 'MPEG-4');
vidfile.FrameRate = 3;
open(vidfile);

ax1 = nexttile;
scatter3(ax1, loc(1, :), loc(2, :), loc(3, :), size, ceil(infections / t), 'filled')
view(ax1, 45 + 90, 45)
colormap cool
axis equal
xlim([0 X])    
ylim([0 Y])
zlim([0 Z])
xlabel('x')
ylabel('y')
zlabel('z')
colorbar('Ticks', [0, 1], 'TickLabels', {'Susceptible', 'Infected'})

ax2 = nexttile;
scatter3(ax2, loc(1, :), loc(2, :), loc(3, :), size, ceil(infections / t), 'filled')
view(ax2, 45, 45)
colormap cool
axis equal
xlim([0 X])    
ylim([0 Y])
zlim([0 Z])
xlabel('x')
ylabel('y')
zlabel('z')
pause
for i = 2 : T
    infections = max(0, infections - 1); %allow infection to pass
    theta = 2 * pi * rand(1, N);
    phi = pi * rand(1, N);
    if reflecting
        loc(1, :) = abs(mod(loc(1, :) + D * cos(theta) .* sin(phi) - X,...
            2 * X) - X);
        loc(2, :) = abs(mod(loc(2, :) + D * sin(theta) .* sin(phi) - Y,...
            2 * Y) - Y);
        loc(3, :) = abs(mod(loc(3, :) + D * cos(phi) - Z,...
            2 * Z) - Z);
    else
        loc(1, :) = mod(loc(1, :) + D * cos(theta) .* sin(phi), X);
        loc(2, :) = mod(loc(2, :) + D * sin(theta) .* sin(phi), Y);
        loc(3, :) = mod(loc(3, :) + D * cos(phi), Z);
    end
    for j = find(infections > 0) % go through infected agents
        distances = sqrt((loc(1, :) - loc(1, j)) .^ 2 + (loc(2, :) -...
            loc(2, j)) .^ 2);
        % check which were infected (but remove already infected)
        infections((distances < R) & (infections == 0)) = t; 
    end
    numberinfected(i) = sum(infections > 0);
    
    tile = tiledlayout(1, 2);
    title(tile, ['Day ', num2str(i), ' (Percent infected: ', num2str(100 * numberinfected(i) / N), '%)'])
    
    ax1 = nexttile;
    scatter3(ax1, loc(1, :), loc(2, :), loc(3, :), size, ceil(infections / t), 'filled')
    view(ax1, 45 + 90, 45)
    axis equal
    xlim([0 X])    
    ylim([0 Y])
    zlim([0 Z])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    colorbar('Ticks', [0, 1], 'TickLabels', {'Susceptible', 'Infected'})
    
    ax2 = nexttile;
    scatter3(ax2, loc(1, :), loc(2, :), loc(3, :), size, ceil(infections / t), 'filled')
    view(ax2, 45, 45)
    axis equal
    xlim([0 X])    
    ylim([0 Y])
    zlim([0 Z])
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
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
