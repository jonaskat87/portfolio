function [ ] = simulation2surfaceplots(Tend, actual)

mq = [0.1 : 0.1 : 1 2 3 4 5 10 15 20];
range = 1 : length(mq);
data = zeros(length(mq), 7, 9);
for i = range
    data(i, :, :) = simulation2extraction(9, i, Tend, actual);
    disp(['Finished h=', num2str(mq(i)), 'm.'])
end
[X, Y] = meshgrid(mq, data(1, 1, :));

figure(1)
tiledlayout(2, 1)
nexttile
grid on
surf(X, Y, squeeze(data(:, 2, :) - data(:, 1, :)).')
set(gca, 'Xscale', 'log', 'Yscale', 'log')
colormap autumn
shading interp
title('Time to solidification across domain')
xlabel('$h\hspace{0.1cm}[m]$', 'interpreter', 'latex')
ylabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
zlabel('$(t_{solidification}-t_{emp})/t_{h}$', 'interpreter', 'latex')

nexttile
grid on
surf(X, Y, squeeze(data(:, 3, :)).')
set(gca, 'Xscale', 'log', 'Yscale', 'log')
colormap autumn
shading interp
title('Time to solidification for first lobe only')
xlabel('$h\hspace{0.1cm}[m]$', 'interpreter', 'latex')
ylabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
zlabel('$t_{solidification}/t_{h}$', 'interpreter', 'latex')

figure(2)
grid on
surf(X, Y, squeeze(data(:, 4, :) - data(:, 1, :)).')
set(gca, 'Xscale', 'log', 'Yscale', 'log')
colormap autumn
shading interp
title(['Time for first lobe to reach ', num2str(int64(Tend - 273.15)), '°C'])
xlabel('$h\hspace{0.1cm}[m]$', 'interpreter', 'latex')
ylabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
zlabel(['$(t_{', num2str(int64(Tend - 273.15)), '^{\circ}C}-t_{emp})/t_{h}$'], 'interpreter', 'latex')

figure(3)
tiledlayout(2, 1)
nexttile
grid on
surf(X, Y, squeeze(data(:, 5, :)).')
set(gca, 'Xscale', 'log', 'Yscale', 'log')
colormap autumn
shading interp
title('Height where solidification completed across domain')
xlabel('$h\hspace{0.1cm}[m]$', 'interpreter', 'latex')
ylabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
zlabel('$z/h$', 'interpreter', 'latex')

nexttile
grid on
surf(X, Y, squeeze(data(:, 6, :)).')
set(gca, 'Xscale', 'log', 'Yscale', 'log')
colormap autumn
shading interp
title('Height where solidification completed in the first lobe')
xlabel('$h\hspace{0.1cm}[m]$', 'interpreter', 'latex')
ylabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
zlabel('$z/h$', 'interpreter', 'latex')

figure(4)
grid on
surf(X, Y, squeeze(data(:, 7, :)).')
set(gca, 'Xscale', 'log', 'Yscale', 'log')
colormap autumn
shading interp
title(['Height where first lobe last reached ', num2str(int64(Tend - 273.15)), '°C'])
xlabel('$h\hspace{0.1cm}[m]$', 'interpreter', 'latex')
ylabel('$t_{emp}/t_{h}$', 'interpreter', 'latex')
zlabel(['$z_{', num2str(int64(Tend - 273.15)), '^{\circ}C}/h$'], 'interpreter', 'latex')