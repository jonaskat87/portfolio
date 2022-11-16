function [] = honsimulation(alpha, outsidetemp, height0, dtq, tend, flag)
% flag: 1 or true for plotting at each timestep saved 
% dtq gives the sampling interval (in seconds)
% height is the size of each lava slab
% tend gives the final time to find the crossing time/position

d = 1;

%parameters
scaling = 4;
height = scaling * height0;
scale = 10;
%force solid at top and bottom boundaries (i.e. phi = 1)
%initial condition: phi = 0 on the inside
T0 = 1142+273.15;
outsidetemp = outsidetemp+273.15; %above flow
groundtemp = T0; %below flow (set to be same temperature as lava)

beta = 5.6e-8; % kinetic coefficient [m Pa / K^2]
Tm = 1070+273.15; % melting temperature [Kelvin]
L = 4e5; % latent heat of fusion [J/m^3]
E = 0.5; % interfacial energy [J/m^2]
cp = 2.57e6; % specific heat at constant pressure [J/m^3*K]
tau = ((d ^ 2) * L / (8 * (log(3) ^ 2) * E * Tm)) * (beta + (209 * d * L / (1680 * log(3) * alpha * cp)));
H = 12 * log(3) * E / d; % energy barrier for solidification
omega = sqrt(2) * d / log(81); % square root of interfacial coefficient [sqrt(m)]
v = 2; % wind speed [m/s]

%scheme-related
dx = 2.5e-3; %meters
dt0 = 1e-11;
dt = dt0;
% m = 5e3; %(make arbitrarily large)
m = 1e4; %(make arbitrarily large)
tq = dtq * (1 : m);
error = 1e-11;
t = 0;
n = round(height / dx);
milne = 270 / 19;

tic
D = sparse(1 : n, 1 : n, -30 * ones(1, n) / 12, n, n);
E = sparse(2 : n, 1 : n - 1, 16 * ones(1, n - 1) / 12, n, n);
F = sparse(3 : n, 1 : n - 2, -ones(1, n - 2) / 12, n, n);
fdstencil = F + E + D + E' + F'; %fourth-order 2nd derivative central difference stencil

fphi = @(phi, T) (((omega / dx) ^ 2) * (fdstencil * phi + [15 / 12; -1 / 12; zeros(n - 4, 1);...
    -1 / 12; 15 / 12]) + 2 * phi .* (1 - phi) .* (2 * phi - 1 - (3 * L / (H * Tm)) * (T - Tm))) / tau;
phi = zeros(n, m);
T = T0 * ones(n, m);
T(1 : int64(n - (n / scaling)), :) = groundtemp * ones(int64(n - (n / scaling)), m);

%phase = zeros(length(xq) * int64(startflow / dt), 2);
const = alpha * cp;
% A = height0 ^ 2;
sigma = 5.67051e-8;
eps = 0.6;
hc = 12.12-1.16*v+11.6*sqrt(v);
% constconv = 60 * dx * A * hc / const;
constconv = 60 * dx * hc / const;
constrad = 60 * sigma * eps * dx / const;
neumann = @(T) constconv * (outsidetemp - T) + constrad * ((outsidetemp ^ 4) - (T ^ 4));
TOne = @(T) (-3 * T(end - 4) + 20 * T(end - 3) - 60 * T(end - 2) + 120 * T(end - 1) - 65 * T(end) + neumann(T(end))) / 12;
TTwo = @(T) (30 * TOne(T) + 20 * T(end) - 60 * T(end - 1) + 15 * T(end - 2) - 2 * T(end - 3) - neumann(T(end))) / 3; 
fT = @(phi, Tl) (alpha / (dx ^ 2)) * (fdstencil * Tl + [15 * groundtemp / 12; -groundtemp / 12;...
    zeros(n - 4, 1); -TOne(Tl) / 12 ; ((16 * TOne(Tl)) - TTwo(Tl)) / 12]) + ((6 * L / cp) * phi .* (1 - phi) .* fphi(phi, Tl));

fphiq = zeros(n, 4);
fTq = zeros(n, 4);

xq = dx * (0 : round((height / dx) / scale) : (height / dx));
startcounter = 1;
i = 1;
k = 1;
Ti = T(:, 1);
phii = phi(:, 1);
firsttime = 1;
while tq(k) < tend
    if (norm(phii - 1, Inf) < error) && (tq(k - 1) == t)
        if firsttime
           phii(:) = 1;
           firsttime = 0;
        end
        if startcounter < 4
            fphiq(:, startcounter) = fphi(phii, Ti);
            fTq(:, startcounter) = fT(phii, Ti);
            k1phi = dt * fphiq(:, startcounter);
            k1T = dt * fTq(:, startcounter);
            k2phi = dt * fphi(phii + 0.4 * k1phi, Ti + 0.4 * k1T);
            k2T = dt * fT(phii + 0.4 * k1phi, Ti + 0.4 * k1T);
            k3phi = dt * fphi(phii + 0.29697761 * k1phi + 0.15875964 * k2phi, Ti + 0.29697761 * k1T + 0.15875964 * k2T);
            k3T = dt * fT(phii + 0.29697761 * k1phi + 0.15875964 * k2phi, Ti + 0.29697761 * k1T + 0.15875964 * k2T);
            k4phi = dt * fphi(phii + 0.2181004 * k1phi - 3.05096516 * k2phi + 3.83286476 * k3phi,...
                Ti + 0.2181004 * k1T - 3.05096516 * k2T + 3.83286476 * k3T);
            k4T = dt * fT(phii + 0.2181004 * k1phi - 3.05096516 * k2phi + 3.83286476 * k3phi,...
                Ti + 0.2181004 * k1T - 3.05096516 * k2T + 3.83286476 * k3T);
            Ti = Ti + 0.17476028 * k1T - 0.55148066 * k2T + 1.2055356 * k3T + 0.17118478 * k4T;
            startcounter = startcounter + 1;    
        else
            fphiq(:, 4) = fphi(phii, Ti);
            fTq(:, 4) = fT(phii, Ti);
            phitilde = phii + (dt / 24) * (55 * fphiq(:, 4) - 59 * fphiq(:, 3) + 37 * fphiq(:, 2) - 9 * fphiq(:, 1));
            Ttilde = Ti + (dt / 24) * (55 * fTq(:, 4) - 59 * fTq(:, 3) + 37 * fTq(:, 2) - 9 * fTq(:, 1));
            phitilde1 = phii + (dt / 24) * (9 * fphi(phitilde, Ttilde) + 19 * fphiq(:, 4) - 5 * fphiq(:, 3) + fphiq(:, 2));
            Ttilde1 = Ti + (dt / 24) * (9 * fT(phitilde, Ttilde) + 19 * fTq(:, 4) - 5 * fTq(:, 3) + fTq(:, 2));
            for iter = 1 : 3
                phitilde1 = phii + (dt / 24) * (9 * fphi(phitilde1, Ttilde1) + 19 * fphiq(:, 4) - 5 * fphiq(:, 3) + fphiq(:, 2));
                Ttilde1 = Ti + (dt / 24) * (9 * fT(phitilde1, Ttilde1) + 19 * fTq(:, 4) - 5 * fTq(:, 3) + fTq(:, 2));
            end
            Tf = Ti + (dt / 24) * (9 * fT(phitilde1, Ttilde1) + 19 * fTq(:, 4) - 5 * fTq(:, 3) + fTq(:, 2));

            %adaptive time stepping
            qT = (error * milne / (sqrt(dx) * norm(Tf - Ttilde))) ^ (1 / 4);

            if qT < 0.8
                dt = dt / 2;
                Ti = Tf;
                startcounter = 1;
                %once stepsize is changed, go back and do RK4 again
            elseif qT > 2
                dt = min(2 * dt, 0.5 * (dx ^ 2) / alpha);
                Ti = Tf;
                startcounter = 1;
                %once stepsize is changed, go back and do RK4 again
            else
                Ti = Tf;
                fphiq = fphiq(:, [2 3 4 1]); 
                fTq = fTq(:, [2 3 4 1]);
            end
        end
        t = t + dt;
        i = i + 1;
        if tq(k) == t
            phi(:, k + 1) = phii;
            T(:, k + 1) = Ti;
            i = 1;
            k = k + 1;
            if flag == 1
                figure(1)
                tiledlayout(2,1)
                % Top plot
                ax1 = nexttile;
                plot(ax1, phii)
                title(ax1, ['Hours elapsed: ' num2str(t / (60 * 60))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                % Bottom plot
                ax2 = nexttile;
                plot(ax2, Ti - 273.15)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [°C]') 
                ylim([0 T0 * 1.1])             
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))          
            end
        elseif t > tq(k)
            t = t - dt;
            dt = tq(k) - t;
            t = t + dt;
            i = 1;
            if startcounter == 1
                fphiq(:, 1) = fphiq(:, 4);
                fTq(:, 1) = fTq(:, 4);
            else
                fphiq(:, 1) = fphiq(:, startcounter - 1);
                fTq(:, 1) = fTq(:, startcounter - 1);        
            end
            k1phi = dt * fphiq(:, 1);
            k1T = dt * fTq(:, 1);
            k2phi = dt * fphi(phii + 0.4 * k1phi, Ti + 0.4 * k1T);
            k2T = dt * fT(phii + 0.4 * k1phi, Ti + 0.4 * k1T);
            k3phi = dt * fphi(phii + 0.29697761 * k1phi + 0.15875964 * k2phi, Ti + 0.29697761 * k1T + 0.15875964 * k2T);
            k3T = dt * fT(phii + 0.29697761 * k1phi + 0.15875964 * k2phi, Ti + 0.29697761 * k1T + 0.15875964 * k2T);
            k4phi = dt * fphi(phii + 0.2181004 * k1phi - 3.05096516 * k2phi + 3.83286476 * k3phi,...
                Ti + 0.2181004 * k1T - 3.05096516 * k2T + 3.83286476 * k3T);
            k4T = dt * fT(phii + 0.2181004 * k1phi - 3.05096516 * k2phi + 3.83286476 * k3phi,...
                Ti + 0.2181004 * k1T - 3.05096516 * k2T + 3.83286476 * k3T);
            Ti = Ti + 0.17476028 * k1T - 0.55148066 * k2T + 1.2055356 * k3T + 0.17118478 * k4T;
            phi(:, k + 1) = phii;
            T(:, k + 1) = Ti;
            startcounter = 2;
            if flag == 1
                figure(1)
                tiledlayout(2,1)
                ax1 = nexttile;
                plot(ax1, phii)
                title(ax1, ['Hours elapsed: ' num2str(t / (60 * 60))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                ax2 = nexttile;
                plot(ax2, Ti - 273.15)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [°C]') 
                ylim([0 T0 * 1.1])             
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))        
            end
            k = k + 1;
        end
    else
        %AB4/AM4 with Ralston's RK4 to initialize
        if startcounter < 4
            fphiq(:, startcounter) = fphi(phii, Ti);
            fTq(:, startcounter) = fT(phii, Ti);
            k1phi = dt * fphiq(:, startcounter);
            k1T = dt * fTq(:, startcounter);
            k2phi = dt * fphi(phii + 0.4 * k1phi, Ti + 0.4 * k1T);
            k2T = dt * fT(phii + 0.4 * k1phi, Ti + 0.4 * k1T);
            k3phi = dt * fphi(phii + 0.29697761 * k1phi + 0.15875964 * k2phi, Ti + 0.29697761 * k1T + 0.15875964 * k2T);
            k3T = dt * fT(phii + 0.29697761 * k1phi + 0.15875964 * k2phi, Ti + 0.29697761 * k1T + 0.15875964 * k2T);
            k4phi = dt * fphi(phii + 0.2181004 * k1phi - 3.05096516 * k2phi + 3.83286476 * k3phi,...
                Ti + 0.2181004 * k1T - 3.05096516 * k2T + 3.83286476 * k3T);
            k4T = dt * fT(phii + 0.2181004 * k1phi - 3.05096516 * k2phi + 3.83286476 * k3phi,...
                Ti + 0.2181004 * k1T - 3.05096516 * k2T + 3.83286476 * k3T);
            phii = phii + 0.17476028 * k1phi - 0.55148066 * k2phi + 1.2055356 * k3phi + 0.17118478 * k4phi;
            Ti = Ti + 0.17476028 * k1T - 0.55148066 * k2T + 1.2055356 * k3T + 0.17118478 * k4T;
            startcounter = startcounter + 1;    
       else
            fphiq(:, 4) = fphi(phii, Ti);
            fTq(:, 4) = fT(phii, Ti);
            phitilde = phii + (dt / 24) * (55 * fphiq(:, 4) - 59 * fphiq(:, 3) + 37 * fphiq(:, 2) - 9 * fphiq(:, 1));
            Ttilde = Ti + (dt / 24) * (55 * fTq(:, 4) - 59 * fTq(:, 3) + 37 * fTq(:, 2) - 9 * fTq(:, 1));
            phitilde1 = phii + (dt / 24) * (9 * fphi(phitilde, Ttilde) + 19 * fphiq(:, 4) - 5 * fphiq(:, 3) + fphiq(:, 2));
            Ttilde1 = Ti + (dt / 24) * (9 * fT(phitilde, Ttilde) + 19 * fTq(:, 4) - 5 * fTq(:, 3) + fTq(:, 2));
            for iter = 1 : 3
                phitilde1 = phii + (dt / 24) * (9 * fphi(phitilde1, Ttilde1) + 19 * fphiq(:, 4) - 5 * fphiq(:, 3) + fphiq(:, 2));
                Ttilde1 = Ti + (dt / 24) * (9 * fT(phitilde1, Ttilde1) + 19 * fTq(:, 4) - 5 * fTq(:, 3) + fTq(:, 2));
            end
            phif = phii + (dt / 24) * (9 * fphi(phitilde1, Ttilde1) + 19 * fphiq(:, 4) - 5 * fphiq(:, 3) + fphiq(:, 2));
            Tf = Ti + (dt / 24) * (9 * fT(phitilde1, Ttilde1) + 19 * fTq(:, 4) - 5 * fTq(:, 3) + fTq(:, 2));

            %adaptive time stepping
            qphi = (error * milne / (sqrt(dx) * norm(phif - phitilde))) ^ (1 / 4);
            qT = (error * milne / (sqrt(dx) * norm(Tf - Ttilde))) ^ (1 / 4);

            if (qphi < 0.8) || (qT < 0.8)
                dt = dt / 2;
                phii = phif;
                Ti = Tf;
                startcounter = 1;
                %once stepsize is changed, go back and do RK4 again
            elseif (qphi > 2) && (qT > 2)
                dt = min(2 * dt, 0.5 * (dx ^ 2) / alpha);
                phii = phif;
                Ti = Tf;
                startcounter = 1;
                %once stepsize is changed, go back and do RK4 again
            else
                phii = phif;
                Ti = Tf;
                fphiq = fphiq(:, [2 3 4 1]); 
                fTq = fTq(:, [2 3 4 1]);
            end
        end
        t = t + dt;
        i = i + 1;
        if tq(k) == t
            phi(:, k + 1) = phii;
            T(:, k + 1) = Ti;
            i = 1;
            k = k + 1;
            if flag == 1
                figure(1)
                tiledlayout(2,1)
                % Top plot
                ax1 = nexttile;
                plot(ax1, phii)
                title(ax1, ['Hours elapsed: ' num2str(t / (60 * 60))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                % Bottom plot
                ax2 = nexttile;
                plot(ax2, Ti - 273.15)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [°C]') 
                ylim([0 T0 * 1.1])             
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))
            else
                disp(t / (60 * 60)) %in hours           
            end
        elseif t > tq(k)
            t = t - dt;
            dt = tq(k) - t;
            t = t + dt;
            i = 1;
            if startcounter == 1
                fphiq(:, 1) = fphiq(:, 4);
                fTq(:, 1) = fTq(:, 4);
            else
                fphiq(:, 1) = fphiq(:, startcounter - 1);
                fTq(:, 1) = fTq(:, startcounter - 1);        
            end
            k1phi = dt * fphiq(:, 1);
            k1T = dt * fTq(:, 1);
            k2phi = dt * fphi(phii + 0.4 * k1phi, Ti + 0.4 * k1T);
            k2T = dt * fT(phii + 0.4 * k1phi, Ti + 0.4 * k1T);
            k3phi = dt * fphi(phii + 0.29697761 * k1phi + 0.15875964 * k2phi, Ti + 0.29697761 * k1T + 0.15875964 * k2T);
            k3T = dt * fT(phii + 0.29697761 * k1phi + 0.15875964 * k2phi, Ti + 0.29697761 * k1T + 0.15875964 * k2T);
            k4phi = dt * fphi(phii + 0.2181004 * k1phi - 3.05096516 * k2phi + 3.83286476 * k3phi,...
                Ti + 0.2181004 * k1T - 3.05096516 * k2T + 3.83286476 * k3T);
            k4T = dt * fT(phii + 0.2181004 * k1phi - 3.05096516 * k2phi + 3.83286476 * k3phi,...
                Ti + 0.2181004 * k1T - 3.05096516 * k2T + 3.83286476 * k3T);
            phii = phii + 0.17476028 * k1phi - 0.55148066 * k2phi + 1.2055356 * k3phi + 0.17118478 * k4phi;
            Ti = Ti + 0.17476028 * k1T - 0.55148066 * k2T + 1.2055356 * k3T + 0.17118478 * k4T;
            phi(:, k + 1) = phii;
            T(:, k + 1) = Ti;
            startcounter = 2;
            if flag == 1
                figure(1)
                tiledlayout(2,1)
                ax1 = nexttile;
                plot(ax1, phii)
                title(ax1, ['Hours elapsed: ' num2str(t / (60 * 60))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                ax2 = nexttile;
                plot(ax2, Ti - 273.15)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [°C]') 
                ylim([0 T0 * 1.1])             
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq')) 
            else
                disp(t / (60 * 60)) %in hours           
            end
            k = k + 1;
        end
        if any(phii > 1.01) || any(phii < -0.01)
            disp('RK4 failed. Decreasing stepsize.')
            t = t - dt;
            i = i - 1;
            dt = dt / 2;
            startcounter = 1;
            phii = phii - 0.17476028 * k1phi + 0.55148066 * k2phi - 1.2055356 * k3phi - 0.17118478 * k4phi;
            Ti = Ti - 0.17476028 * k1T + 0.55148066 * k2T - 1.2055356 * k3T - 0.17118478 * k4T;
        end
    end
end
toc;
disp(['The MOL part of the simulation completed in ',...
    num2str(toc / 60), ' minutes.'])

% phi = phi(:, 1 : k);
T = T(:, 1 : k);
tq = tq(1 : k);

T800 = zeros(1, k);
T1070 = zeros(1, k);
dist800 = abs(T - (800 + 273.15));
dist1070 = abs(T - (1070 + 273.15));
xqfull = (1 : n) * dx;
for i = 2 : k
    [~, I800] = min(dist800(:, i));
    [~, I1070] = min(dist1070(:, i));
    T800(i) = height - xqfull(I800);
    T1070(i) = height - xqfull(I1070);
end

% figure(3)
% tiledlayout(2,1)
% ax3 = nexttile;
% surf(ax3, phi)
% colormap autumn
% shading interp
% xticks(0 : k / scale : k)
% xticklabels(num2str(round(t / (60 * 60 * k), 1) * (0 : k / scale : k)'))
% xlabel('time [hours]')
% yticks(0 : n / scale : n)
% yticklabels(num2str(xq'))
% ylabel('height [m]')
% zlabel('Φ (0 for liquid, 1 for solid)')
% 
% ax4 = nexttile;
% surf(ax4, T - 273.15)
% colormap autumn
% shading interp
% xticks(0 : k / scale : k)
% xticklabels(num2str(round(t / (60 * 60 * k), 1) * (0 : k / scale : k)'))
% xlabel('time [hours]')
% yticks(0 : n / scale : n)
% yticklabels(num2str(xq'))
% ylabel('height [m]')
% zlabel('T [°C]')

% cc = flipud(hot(k));
% for index = 1 : k
%     figure(4); hold on;
%     plot(T(:, index), phi(:, index), 'o', 'color', cc(index, :));
% end

load hondata.mat fig9 fig10  
data800 = fig9.colder;
data1070 = fig9.hotter;
s = fig10.s;
a = fig10.a;
b = fig10.b;
c = fig10.c;
d = fig10.d;
e = fig10.e;

colores = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560; 0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330];

figure(5)
% crust = mean(phi(int64((height - height0) / dx) : end, :), 1);
% figure(5); plot(sqrt(tq / (60 * 60)), crust, 'ro');
h = plot(sqrt(tq(2 : end) / (60 * 60)), [T800(2 : end); T1070(2 : end)]);
for i = 1 : 2
    set(h(i), 'color', colores(i, :)) 
end
hold on
scatter(data800(:, 1), data800(:, 2), [], colores(1, :))
scatter(data1070(:, 1), data1070(:, 2), [], colores(2, :))
hold off
legend({'800°C', '1070°C'}, 'Location', 'northwest')
xlabel('sqrt(time) [sqrt(hours)]')
ylabel('Thickness of crust [m]')
set(gca, 'fontsize', 14)

tempsample = n - [0 round(0.005 / dx) round(0.02 / dx) round(0.06 / dx) round(0.1 / dx) round(0.2 / dx)];

figure(6)
h = plot(tq(2 : end) / (60 * 60), T(tempsample, 2 : end) - 273.15);
for i = 1 : length(tempsample)
   set(h(i), 'color', colores(i, :)) 
end
hold on
scatter(s(:, 1), s(:, 2), [], colores(1, :))
scatter(a(:, 1), a(:, 2), [], colores(2, :))
scatter(b(:, 1), b(:, 2), [], colores(3, :))
scatter(c(:, 1), c(:, 2), [], colores(4, :))
scatter(d(:, 1), d(:, 2), [], colores(5, :))
scatter(e(:, 1), e(:, 2), [], colores(6, :))
hold off
legend({'surface', '0.005 m', '0.02 m', '0.06 m', '0.10 m', '0.20 m'}, 'Location', 'northeast')
xlabel('time [hours]')
ylabel('T [°C]')
set(gca, 'fontsize', 14)

figure(7)
h = semilogx(tq(2 : end) / (60 * 60), T(tempsample, 2 : end) - 273.15);
for i = 1 : length(tempsample)
   set(h(i), 'color', colores(i, :)) 
end
hold on
scatter(s(:, 1), s(:, 2), [], colores(1, :))
scatter(a(:, 1), a(:, 2), [], colores(2, :))
scatter(b(:, 1), b(:, 2), [], colores(3, :))
scatter(c(:, 1), c(:, 2), [], colores(4, :))
scatter(d(:, 1), d(:, 2), [], colores(5, :))
scatter(e(:, 1), e(:, 2), [], colores(6, :))
hold off
legend({'surface', '0.005 m', '0.02 m', '0.06 m', '0.10 m', '0.20 m'}, 'Location', 'southwest')
xlabel('time [hours]')
ylabel('T [°C]')
set(gca, 'fontsize', 14)