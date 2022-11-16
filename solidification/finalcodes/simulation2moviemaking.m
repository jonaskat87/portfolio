function [ ] = simulation2moviemaking(i, height0, height1, Tend)
cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\newemplacementresults'

Tm = 1070+273.15;
T0 = 1200+273.15;

% change this manually
fname = ['emplacementresults_', num2str(height0), '_', num2str(height1),...
 '_', num2str(round(Tend)), 'K.mat'];
m = matfile(fname);

all = m.(['all', num2str(i)]); 
phi = all{2};
T = all{3};
dx = all{4};
k = all{5};
times = all{7};
t = times(2);
n = length(T(:, 1));
tq = (t / k) * (1 : k);

moviename = ['emplacementresults_', num2str(height0), '_', num2str(height1),...
 '_', num2str(round(Tend)), 'K_',...
 num2str(round(times(1) / (60 * 60))), 'hours.mp4'];

fh = figure(1);
fh.WindowState = 'maximized';
vidfile = VideoWriter(moviename, 'MPEG-4');
vidfile.FrameRate = 10;
open(vidfile);
for i = 1 : k
    tiledlayout(1, 2)
    ax1 = nexttile; 
    if tq(i) > times(1)
        plot(ax1, phi(:, i), dx * (1 : n) - 4 * height0)
    else
        plot(ax1, phi(1 : (5 * end / 6), i), dx * (1 : (5 * n / 6)) - 4 * height0)    
    end
    hold on
    x1a = xline(1, '--k', 'Solid');
    x1a.LabelVerticalAlignment = 'bottom';
    x1a.LabelOrientation = 'horizontal';
    x1b = xline(0, '--k', 'Liquid');
    x1b.LabelVerticalAlignment = 'bottom';
    x1b.LabelOrientation = 'horizontal';
    hold off
    grid on
    ylabel(ax1, 'Height [m]')
    xlabel(ax1, 'Φ (0 for liquid, 1 for solid)') 
    title(ax1, ['Days elapsed: ' num2str(round(tq(i) / (60 * 60 * 24)))])
    xlim([-0.1 1 + 0.1]) 
    xticks(0 : 0.1 : 1)      
    
    ax2 = nexttile;
    if tq(i) > times(1)
        plot(ax2, T(:, i) - 273.15, dx * (1 : n) - 4 * height0)
    else
        plot(ax2, T(1 : (5 * end / 6), i) - 273.15, dx * (1 : (5 * n / 6)) - 4 * height0)    
    end
    hold on
    x2 = xline(T0 - 273.15, '--k', '$T_{0}$', 'interpreter', 'latex');
    x2.LabelVerticalAlignment = 'bottom';
    x2.LabelOrientation = 'horizontal';
    x3 = xline(Tm - 273.15, '--k', '$T_{m}$', 'interpreter', 'latex');
    x3.LabelVerticalAlignment = 'bottom';
    x3.LabelOrientation = 'horizontal';
    hold off
    grid on   
    ylabel(ax2, 'Height [m]')            
    xlabel(ax2, 'Temperature [°C]') 
    xlim([0 round((T0 - 273.15) * 1.1)])
    xticks((0 : 1 : 8) * round((T0 - 273.15) / 8)) 
    
    drawnow
    F = getframe(gcf); 
    writeVideo(vidfile, F);
end
close(vidfile)