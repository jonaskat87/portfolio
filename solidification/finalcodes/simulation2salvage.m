function [ ] = simulation2salvage(N, height0, height1, Tend)
% loads data and finds temperature checkpoints, as long as the desired
% temperature is BEFORE the one in the data. Creates new file with these.

cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\newemplacementresults'
fname = ['emplacementresults_', num2str(height0), '_', num2str(height1),...
 '_', num2str(round(Tend)), 'K.mat'];
m = matfile(fname, 'Writable', true);
Tm = 1070+273.15;

for i = 1 : N
    all = m.(['all', num2str(i)]);  
    dx = all{4};
    k = all{5};
    heights = all{6};
    height = heights(1);
    phi = all{2};
    T = all{3};
    Tselect = T(int64((height - height0) / dx) + 1 : int64(height / dx), :);
    phiselect = phi(int64((height - height0) / dx) + 1 : int64(height / dx), :);
    times = all{7};
    t = times(2);
    tq = (t / k) * (0 : k - 1);
    xvec = dx * (1 : ((height + heights(3)) / dx));
    xsmallvec = dx * (1 : (height0 / dx));
    flag1 = 0;
    flag2 = 0;
    flag3 = 0;
    for j = flip(1 : k)
        if ~((sum(Tselect(:, j) > Tend) == 0) || flag1)
            tneeded = mean([tq(j) tq(min(k, j + 1))]);
            [~, I] = max(Tselect(:, j));
            xneeded = xsmallvec(min(length(xsmallvec), I));
            all{9} = [tneeded xneeded];
            m.(['all', num2str(i)]) = all;
            flag1 = 1;
        end
        % check for which point was the last to solidify
        if ~((sum(T(:, j) > Tm) == 0) || flag2)
            tlastoverall =  mean([tq(j) tq(min(k, j + 1))]);
            [~, I] = min(phi(:, j));
            xlastoverall = xvec(I) - xvec(int64((height - height0) / dx));
            flag2 = 1;
            all{11} = [tlastoverall xlastoverall];
            m.(['all', num2str(i)]) = all;
        end
        if ~((sum(Tselect(:, j) > Tm) == 0) || flag3)
            tlastfirst =  mean([tq(j) tq(min(k, j + 1))]);
            [~, I] = min(phiselect(:, j));
            xlastfirst = xsmallvec(min(length(xsmallvec), I));
            flag3 = 1;
            all{10} = [tlastfirst xlastfirst];
            m.(['all', num2str(i)]) = all;
        end
        if flag1 && flag2 && flag3 
            break
        end
    end
end

cd ..