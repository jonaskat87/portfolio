mq = [0.1 : 0.1 : 1 2 3 4 5 10 15 20];
ratios = 2 .^ (-4 : 4);
for j = 1 : length(mq)
    h = mq(j);
    d = 1;
    height0 = h;
    height1 = h;
    scaling = 5;
    dx = min(round(sqrt(10) ^ (log10(h) - 5), 3), 0.01); %meters
    Tend = 570 + 273.15;
    dtq = 60 * (h / 0.1) ^ 2;
    
    fname = 'propconstant.mat';
    m = matfile(fname);
    all = m.(['case', num2str(j)]);
    datalastoverall = all{11};
    th = mean(datalastoverall(1, :));
    
    t0 = round(th * ratios(1));
    all1 = simulation2fortestingtemp(scaling, dx, d, height0, height1, t0,...
     dtq, Tend, false);
    fname = ['emplacementresults_', num2str(height0), '_', num2str(height1),...
     '_', num2str(round(Tend)), 'K.mat'];
    save(fname, 'all1');
    m = matfile(fname, 'Writable', true);
    disp(['I just finished the case for a secondary emplacement after ',...
     num2str(t0 / (60 * 60)), ' hours for a ', num2str(h), ' meter lobe.'])

    % run rest of cases
    for i = 2 : length(ratios)
        t = round(th * ratios(i));
        all = simulation2fortestingtemp(scaling, dx, d, height0, height1,...
     t, dtq, Tend, false);
        m.(['all', num2str(i)]) = all;
        disp(['I just finished the case for a secondary emplacement after ',...
     num2str(t / (60 * 60)), ' hours for a ', num2str(h), ' meter lobe.'])
    end
end