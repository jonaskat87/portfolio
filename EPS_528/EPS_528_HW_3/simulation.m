function [ fpclassify ] = simulation(tol, V0, n0, tspan, Vq,...
    I, C, gL, EL, gNa, ENa, Vm, km, gK, EK, Vn, kn, tauV)

lesstol = sqrt(tol);
minf = @(V) 1 ./ (1 + exp((Vm - V) / km));
ninf = @(V) 1 ./ (1 + exp((Vn - V) / kn));
Vfunc = @(V, n) (I - gL * (V - EL) - gNa *...
    minf(V) .* (V - ENa) - gK * n .* (V - EK)) / C;
nfunc = @(V, n) (ninf(V) - n) / tauV;
fun = @(vec) [Vfunc(vec(1), vec(2)); nfunc(vec(1), vec(2))];
Vnullcline = @(V) (I - gL * (V - EL) - gNa *...
    minf(V) .* (V - ENa)) ./ (gK * (V - EK));

% find fixed points by brute force
Vinit = -100 : 1 : 100;
ninit = 0 : 0.1 : 1;
fps = NaN(2, length(Vinit) * length(ninit));
options = optimoptions('fsolve', 'FunctionTolerance', tol,...
    'OptimalityTolerance', tol, 'MaxFunctionEvaluations',...
    1e3, 'Display', 'off');
for i = 1 : length(Vinit)
    for j = 1 : length(ninit)
        % solve for steady states (fixed points)
        trial = fsolve(fun, [Vinit(i) ninit(j)], options);
        % check if actually converged to fixed point
        if norm(fun(trial)) < lesstol 
            fps(:, i + length(Vinit) * (j - 1)) = trial;
        end
    end
end
% remove redundant fixed point values
fps = unique(round(fps, round(-log10(lesstol))).', 'rows').'; 
fps = fps(:, ~(isnan(fps(1, :)) | isnan(fps(2, :))));
disp(fps) % display fixed points located

Jacobian = zeros(2, 2);
% first row is determinant
% second row is trace
fpdata = zeros(2, length(fps(1, :)));
% returns type of fixed point
fpclassify = struct([]);
for i = 1 : length(fps(1, :))
    Vi = fps(1, i);
    ni = fps(2, i);
    % second-order central difference approximation
    % to the first-derivatives
    Jacobian(:, 1) = (fun([Vi + tol ni])...
        - fun([Vi - tol ni])) / (2 * tol);
    Jacobian(:, 2) = (fun([Vi ni + tol])...
        - fun([Vi ni - tol])) / (2 * tol);
    fpdata(1, i) = det(Jacobian);
    fpdata(2, i) = trace(Jacobian);
    % classify fixed point via Jacobian and combinations
    % of determinant and trace values
    if abs(fpdata(1, i)) < lesstol
        if fpdata(2, i) > lesstol
            fpclassify{i} = 'near unstable line';
        elseif fpdata(2, i) < -lesstol
            fpclassify{i} = 'near stable line';            
        else
            fpclassify{i} = 'near uniform motion';
        end   
    elseif (abs(fpdata(2, i)) < lesstol) &&...
            (fpdata(1, i) > lesstol)
        fpclassify{i} = 'near center';
    else
        if fpdata(1, i) < -lesstol
            fpclassify{i} = 'saddle';
        else
            if abs(fpdata(1, i) - 0.25 *...
                    (fpdata(2, i)) ^ 2) < lesstol
                if fpdata(2, i) > lesstol
                    fpclassify{i} = 'near degenerate source';
                else
                    fpclassify{i} = 'near degenerate sink';
                end
            else
                if fpdata(2, i) > lesstol
                    if fpdata(1, i) < 0.25 * (fpdata(2, i)) ^ 2
                        fpclassify{i} = 'source';
                    else
                        fpclassify{i} = 'spiral source';
                    end
                else
                    if fpdata(1, i) < 0.25 * (fpdata(2, i)) ^ 2
                        fpclassify{i} = 'sink';
                    else
                        fpclassify{i} = 'spiral sink';
                    end                    
                end
            end
        end
    end
end

% identify stable, unstable, and other points
stable = (fpdata(1, :) >= 0) & (fpdata(2, :) < -lesstol);
degenerate = (fpdata(1, :) >= 0) & (abs(fpdata(2, :)) < lesstol);
unstable = setdiff(1 : length(fps(1, :)), stable | degenerate);

figure
% plot nullclines
h = zeros(1, 2);
h(1) = plot(Vq, Vnullcline(Vq), 'c--',...
    'DisplayName', 'V-nullcline');
hold on
grid on
h(2) = plot(Vq, ninf(Vq), 'm--',...
    'DisplayName', 'n-nullcline');
title(['Phase portrait for neuron model, I=', num2str(I)],...
    'fontsize', 16)
xlabel('$V(t)$', 'interpreter', 'latex',...
    'fontsize', 16)
ylabel('$n(t)$', 'interpreter', 'latex',...
    'fontsize', 16)
xlim([Vq(1) Vq(end)])
ylim([0 1])
% plot trajectories
opts = odeset('RelTol', tol); % set error tolerance
func = @(t, vec) fun(vec);
colors = {[0 0.4470 0.7410], [0.8500 0.3250 0.0980],...
    [0.9290 0.6940 0.1250], [0.4940 0.1840 0.5560], ...
    [0.4660 0.6740 0.1880], [0.3010 0.7450 0.9330],...
    [0.6350 0.0780 0.1840]};
for i = 1 : length(V0)
    u0 = [V0(i) n0(i)]; % set initial values for x_{1} and x_{2}
    sol = ode45(func, tspan, u0, opts); % use ode45 to solve ODE
    u = sol.y;
    plot(u(1, :), u(2, :), 'Color', colors{i})
    arrowh(u(1, :), u(2, :), colors{i}, [],...
        [0.5 : 0.5 : 2.5 3 : 9] .^ 2);     
end
% plot fixed points 
% stable is blue, unstable is red, and all others are green
scatter(fps(1, stable), fps(2, stable), 'filled', 'b')
scatter(fps(1, degenerate), fps(2, degenerate), 'filled', 'g')
scatter(fps(1, unstable), fps(2, unstable), 'filled', 'r')
hold off
legend(h) % add legend only for nullclines