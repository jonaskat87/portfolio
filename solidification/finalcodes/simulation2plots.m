function [ ] = simulation2plots2(N, Tend, actual)

mq = [0.1 : 0.1 : 1 2 3 4 5 10 15 20];
dataq = ([]);
% range = 1 : 7;
% range = 8 : 13;
% range = 14 : length(mq);
range = 1 : 17;
%range = [1 2 7];

firstend = [3 3 3 3 2 2 2 2 3 3 2 2 2 2 2 2 2];
secondend = [5 4 5 5 4 4 4 4 5 5 4 4 4 4 4 4 4];

load('batlow.mat')
size_step = floor(256/length(range));
all_colors = batlow(1:size_step:256,:); %jet(length(range));
for i = range
    dataq{i} = simulation2extraction(N, i, Tend, actual);
    disp(['Finished h=', num2str(mq(i)), 'm.'])
end

figure(1)
Legend = ([]);
k = 1;
for i = range
    data = dataq{i};
    semilogx(data(1, 1 : firstend(i)), data(2, 1 : firstend(i)),...
        '*', 'Color', all_colors(k, :), 'HandleVisibility', 'off')
    hold on
    semilogx(data(1, (firstend(i) + 1) : secondend(i)),...
        data(2, (firstend(i) + 1) : secondend(i)), 'x','Color',...
        all_colors(k, :),...
        'HandleVisibility', 'off')
    semilogx(data(1, secondend(i) + 1 : end),...
        data(2, secondend(i) + 1 : end), '+', 'Color',...
        all_colors(k, :), 'HandleVisibility', 'off') 
    semilogx(nan, nan, 'o', 'Color', all_colors(k, :));
    Legend{k} = ['h=', num2str(mq(i))];
    k = k + 1;
end
hold off
grid on
title('Time to solidification across domain (w/ emplacement time)', 'fontsize', 16)
xlabel('$t_{emp}/t_{h}$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$t_{solidification}/t_{h}$', 'interpreter', 'latex', 'fontsize', 14)
legend(Legend, 'NumColumns',2)

figure(2)
Legend = ([]);
k = 1;
for i = range
    data = dataq{i};
    semilogx(data(1, 1 : firstend(i)), data(2, 1 : firstend(i))...
        - data(1, 1 : firstend(i)), '*','Color',...
        all_colors(k, :), 'HandleVisibility', 'off')
    hold on
    semilogx(data(1, (firstend(i) + 1) : secondend(i)),...
        data(2, (firstend(i) + 1) : secondend(i)) - data(1,...
        (firstend(i) + 1) : secondend(i)), 'x','Color',...
        all_colors(k, :),...
        'HandleVisibility', 'off')
    semilogx(data(1, secondend(i) + 1 : end),...
        data(2, secondend(i) + 1 : end) - ...
        data(1, secondend(i) + 1 : end), '+','Color',...
        all_colors(k, :), 'HandleVisibility', 'off') 
    semilogx(nan, nan, 'o', 'Color', all_colors(k, :));
    Legend{k} = ['h=', num2str(mq(i))];
    k = k + 1;
end
hold off
grid on
title('Time to solidification across domain (w/o emplacement time)', 'fontsize', 16)
xlabel('$t_{emp}/t_{h}$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$(t_{solidification}-t_{emp})/t_{h}$', 'interpreter', 'latex', 'fontsize', 14)
legend(Legend, 'NumColumns',2)

figure(3)
Legend = ([]);
k = 1;
for i = range
    data = dataq{i};
    semilogx(data(1, 1 : firstend(i)), data(3, 1 : firstend(i)),...
        '*','Color', all_colors(k, :), 'HandleVisibility', 'off')
    hold on
    semilogx(data(1, (firstend(i) + 1) : secondend(i)),...
        data(3, (firstend(i) + 1) : secondend(i)), 'x','Color',...
        all_colors(k, :), 'HandleVisibility', 'off')
    semilogx(data(1, secondend(i) + 1 : end),...
        data(3, secondend(i) + 1 : end), '+','Color',...
        all_colors(k, :), 'HandleVisibility', 'off') 
    semilogx(nan, nan, 'o', 'Color', all_colors(k, :));
    Legend{k} = ['h=', num2str(mq(i))];
    k = k + 1;
end
hold off
grid on
title('Time to solidification for first lobe only', 'fontsize', 16)
xlabel('$t_{emp}/t_{h}$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$t_{solidification}/t_{h}$', 'interpreter', 'latex', 'fontsize', 14)
legend(Legend, 'NumColumns',2)

figure(4)
Legend = ([]);
k = 1;
for i = range
    data = dataq{i};
    semilogx(data(1, 1 : firstend(i)), data(5, 1 : firstend(i)),...
        '*','Color', all_colors(k, :), 'HandleVisibility', 'off')
    hold on
    semilogx(data(1, (firstend(i) + 1) : secondend(i)),...
        data(5, (firstend(i) + 1) : secondend(i)), 'x','Color',...
        all_colors(k, :), 'HandleVisibility', 'off')
    semilogx(data(1, secondend(i) + 1 : end),...
        data(5, secondend(i) + 1 : end), '+','Color',...
        all_colors(k, :), 'HandleVisibility', 'off') 
    semilogx(nan, nan, 'o', 'Color', all_colors(k, :));
    Legend{k} = ['h=', num2str(mq(i))];
    k = k + 1;
end
hold off
grid on
title('Height where solidification completed across domain', 'fontsize', 16)
xlabel('$t_{emp}/t_{h}$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$z/h$', 'interpreter', 'latex', 'fontsize', 14)
legend(Legend, 'Location', 'southeast', 'NumColumns',2)

figure(5)
Legend = ([]);
k = 1;
for i = range
    data = dataq{i};
    semilogx(data(1, 1 : firstend(i)), data(6, 1 : firstend(i)),...
        '*','Color', all_colors(k, :), 'HandleVisibility', 'off')
    hold on
    semilogx(data(1, (firstend(i) + 1) : secondend(i)),...
        data(6, (firstend(i) + 1) : secondend(i)), 'x','Color',...
        all_colors(k, :), 'HandleVisibility', 'off')
    semilogx(data(1, secondend(i) + 1 : end),...
        data(6, secondend(i) + 1 : end), '+','Color',...
        all_colors(k, :), 'HandleVisibility', 'off') 
    semilogx(nan, nan, 'o', 'Color', all_colors(k, :));
    Legend{k} = ['h=', num2str(mq(i))];
    k = k + 1;
end
hold off
grid on
title('Height where solidification completed in the first lobe', 'fontsize', 16)
xlabel('$t_{emp}/t_{h}$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$z/h$', 'interpreter', 'latex', 'fontsize', 14)
legend(Legend, 'Location', 'east', 'NumColumns',2)