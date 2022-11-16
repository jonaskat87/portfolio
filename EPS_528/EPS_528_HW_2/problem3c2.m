fps = [-1 0 1; 0 0 0]; % fixed points
spacing = 0.2; % spacing for vector plot
[X1, X2] = meshgrid(-2: spacing: 2, -2: spacing: 2); 
% vector positions
U = X2; % \dot{x_{1}}
V = X1 .* (1 - X1 .^ 2); % \dot{x_{2}}

% sample trajectories
x10q = [-0.8 : 0.2 : 0.8 1.5 : 0.1 : 2]; 
% we take x_{2}(0)=0 for all x_{2} 
t0 = 0; % initial time
tf = 200; % final time

figure(1)
quiver(X1, X2, U, V, 'b') % plot flow field
title('Phase portrait for the Duffing oscillator',...
    'fontsize', 16)
xlabel('$x_{1}=x$', 'interpreter', 'latex',...
    'fontsize', 16)
ylabel('$x_{2}=\dot{x}$', 'interpreter', 'latex',...
    'fontsize', 16)
xlim([-2.1 2.1])
ylim([-2.2 2.2])
hold on
% plot trajectories
for i = 1 : length(x10q)
    sol = problem3c1(x10q(i), 0, t0, tf);
    u = sol.y;
    plot(u(1, :), u(2, :), 'k-')
end
% plot eigenvectors
annotation('line', [0.46 0.5] + 0.02,...
    [0.54 0.5] + 0.02, 'color', 'r')
annotation('line', [0.54 0.5] + 0.02,...
    [0.46 0.5] + 0.02, 'color', 'r')
annotation('arrow', [0.47 0.471] + 0.02,...
    [0.53 0.529] + 0.02, 'color', 'r')
annotation('arrow', [0.53 0.529] + 0.02,...
    [0.47 0.471] + 0.02, 'color', 'r')
annotation('arrow', [0.5 0.55] + 0.02,...
    [0.5 0.55] + 0.02, 'color', 'r')
annotation('arrow', [0.5 0.45] + 0.02,...
    [0.5 0.45] + 0.02, 'color', 'r')
% plot fixed points
scatter(fps(1, :), fps(2, :), 'filled', 'g')
hold off