function [ data ] = BTWfirst(L, nmax, plotting, disp)
% Model size is LxL
% nmax is the number of iterations
% plotting = 1 (or true) if plot intermediate results
% plotting = 0 (or false) if do not plot intermediate results
% set disp to the number of plots you want to display
% from intermediate iterations. Samples these at equally-
% spaced intervals, and does NOT include the initial one.
% if disp = 0, then no intermediate steps are saved.

% iterates to save
iter = 1;
if disp > 0
   iterq = [(1 : disp) * nmax / (disp + 1) nmax + 1];
else
   iterq = nmax + 1;
end

% initialize grid
U = randi(4, L) - 1;
% store avalanche count here
avalanche = zeros(L);
avalanchesizes = zeros(1, nmax);
% iterate simulation
if plotting || (disp > 0)
    figure(1)
    imagesc(U)
    colormap(flipud(autumn))
    colorbar
    caxis([0 4])
    title('Iterations: 0')
    xlabel('x-position')
    ylabel('y-position')
end
for i = 1 : nmax
    ind = randi(L, 1, 2);
    U(ind(1), ind(2)) = U(ind(1), ind(2)) + 1;
    if plotting
        figure(2)
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
        	figure(2)
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
    if len > 0
        randomv = randi(len);
        I = I(randomv);
        J = J(randomv);
    end
    while len > 0
        if I ~= L
            U(I + 1, J) = U(I + 1, J) + 1;
        end
        if I ~= 1
            U(I - 1, J) = U(I - 1, J) + 1;
        end
        if J ~= L
            U(I, J + 1) = U(I, J + 1) + 1;
        end
        if J ~= 1
            U(I, J - 1) = U(I, J - 1) + 1;
        end
        U(I, J) = 0;
        avalanche(I, J) = 1;
        if plotting
            figure(2)
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
        if len > 0
            randomv = randi(len);
            I = I(randomv);
            J = J(randomv);
        end
    end
    if (i >= iterq(iter)) && (disp > 0)
        figure(3 + iter - 1)
        imagesc(U)
        colormap(flipud(autumn))
        colorbar
        caxis([0 4])
        title(['Iterations: ', num2str(i)])
        xlabel('x-position')
        ylabel('y-position')  
        iter = iter + 1;
    end
    avalanchesizes(i) = sum(sum(avalanche));
    avalanche = zeros(L);
end
if plotting || (disp > 0)
    figure(3 + disp)
    imagesc(U)
    colormap(flipud(autumn))
    colorbar
    caxis([0 4])
    title(['Iterations: ', num2str(i)])
    xlabel('x-position')
    ylabel('y-position')
end
[count, domain] = histcounts(avalanchesizes,...
    'BinMethod', 'integers');
domain = domain(1 : end - 1) + 0.5;
data = struct('sizes', avalanchesizes,...
    'occurrence', [count; domain]);