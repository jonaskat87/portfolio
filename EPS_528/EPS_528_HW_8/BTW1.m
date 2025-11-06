function [ data ] = BTWfirst(L, nmax, plotting)
% Model size is LxL
% nmax is the number of iterations
% plotting = 1 (or true) if plot intermediate results
% plotting = 0 (or false) if do not plot intermediate results

% initialize grid
U = randi(4, L) - 1;
% store avalanche count here
avalanche = zeros(L);
avalanchesizes = zeros(1, nmax);
% iterate simulation
i = 0;
if plotting
    figure(1)
    imagesc(U)
    colormap(flipud(autumn))
    colorbar
    caxis([0 4])
    title(['Iterations: ', num2str(i)])
    xlabel('x-position')
    ylabel('y-position')
end
for i = 1 : nmax
    ind = randi(L, 1, 2);
    U(ind(1), ind(2)) = U(ind(1), ind(2)) + 1;
    if plotting
        figure(1)
        imagesc(U)
        colormap(flipud(autumn))
        colorbar
        caxis([0 4])
        title(['Iterations: ', num2str(i)])
        xlabel('x-position')
        ylabel('y-position')    
    end
    % relax first site if needed because we know where it is
    if U(ind(1), ind(2)) > 3
        if ind(1) ~= L
            U(ind(1) + 1, ind(2)) = U(ind(1) + 1, ind(2)) + 1;
        end
        if ind(1) ~= 1
            U(ind(1) - 1, ind(2)) = U(ind(1) - 1, ind(2)) + 1;
        end
        if ind(2) ~= L
            U(ind(1), ind(2) + 1) = U(ind(1), ind(2) + 1) + 1;
        end
        if ind(2) ~= 1
            U(ind(1), ind(2) - 1) = U(ind(1), ind(2) - 1) + 1;
        end
        U(ind(1), ind(2)) = 0;
        avalanche(ind(1), ind(2)) = 1;
        if plotting
        	figure(1)
            imagesc(U)
            colormap(flipud(autumn))
            colorbar
            caxis([0 4])
            title(['Iterations: ', num2str(i)])
            xlabel('x-position')
            ylabel('y-position')  
        end
    end
    % find other unstable locations
    [I, J] = find(U > 3);
    len = length(I);
    while len > 0
        for k = 1 : len
            if I(k) ~= L
                U(I(k) + 1, J(k)) = U(I(k) + 1, J(k)) + 1;
            end
            if I(k) ~= 1
                U(I(k) - 1, J(k)) = U(I(k) - 1, J(k)) + 1;
            end
            if J(k) ~= L
                U(I(k), J(k) + 1) = U(I(k), J(k) + 1) + 1;
            end
            if J(k) ~= 1
                U(I(k), J(k) - 1) = U(I(k), J(k) - 1) + 1;
            end
            U(I(k), J(k)) = 0;
            avalanche(I(k), J(k)) = 1;
        end
        if plotting
            figure(1)
            imagesc(U)
            colormap(flipud(autumn))
            colorbar
            caxis([0 4])
            title(['Iterations: ', num2str(i)])
            xlabel('x-position')
            ylabel('y-position')  
        end
        [I, J] = find(U > 3);
        len = length(I);
    end
    avalanchesizes(i) = sum(sum(avalanche));
    avalanche = zeros(L);
end
[count, domain] = histcounts(avalanchesizes, 'BinMethod', 'integers');
domain = int64(domain(1 : end - 1) + 0.5);
data = struct('sizes', avalanchesizes, 'occurrence', [count; domain]);