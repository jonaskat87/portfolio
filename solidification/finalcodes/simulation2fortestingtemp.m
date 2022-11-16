function [ all ] = simulation2fortestingtemp(scaling, dx, d, height0, height1, t1, dtq, Tend, flag)
% flag: 1 or true for plotting at each timestep saved 
% dtq gives the sampling interval
% height is the size of each lava slab
% Tend gives the final temperature to find the crossing time/position

%parameters
height = scaling * height0;
scale = 10;
%force solid at top and bottom boundaries (i.e. phi = 1)
%initial condition: phi = 0 on the inside
%set outside temperature to 30 C, ground below lava to 20 C
%lava 1600 C
outsidetemp = 30+273.15; %above flow
groundtemp = 20+273.15; %below flow

beta = 5.6e-8; % kinetic coefficient [m Pa / K^2]
Tm = 1070+273.15; % melting temperature [Kelvin]
L = 4e5; % latent heat of fusion [J/m^3]
E = 0.5; % interfacial energy [J/m^2]
alpha = 3.75e-7;  % thermal diffusivity [m^2/s]
cp = 2.57e6; % specific heat at constant pressure [J/m^3*K]
tau = ((d ^ 2) * L / (8 * (log(3) ^ 2) * E * Tm)) * (beta + (209 * d * L / (1680 * log(3) * alpha * cp)));
H = 12 * log(3) * E / d; % energy barrier for solidification
omega = sqrt(2) * d / log(81); % square root of interfacial coefficient [sqrt(m)]
v=2; % wind speed

% L=350; % latent heat of fusion [J/kg]
% H = L/8; % energy barrier for solidification
% cp=800; % specific heat at constant pressure [J/kg*K]
% alpha=3.75e-7;  % thermal diffusivity [m^2/s]
% omega=1e-1; % square root of interfacial coefficient [sqrt(m)]
% Tm=1070+273.15; % melting temperature [Kelvin]
% tau=150*16; % char. time of solidification [s]
% v=2; % wind speed

%scheme-related
dt0 = 1e-11;
dt = dt0;
m = 1e3; %(make arbitrarily large)
tq = dtq * (1 : m);
error = 1e-11;
t = 0;
n = int64(height / dx);
milne = 270 / 19;

tic
D = sparse(1 : n, 1 : n, -30 * ones(1, n) / 12, n, n);
E = sparse(2 : n, 1 : n - 1, 16 * ones(1, n - 1) / 12, n, n);
F = sparse(3 : n, 1 : n - 2, -ones(1, n - 2) / 12, n, n);
fdstencil = F + E + D + E' + F'; %fourth-order 2nd derivative central difference stencil

fphi = @(phi, T) (((omega / dx) ^ 2) * (fdstencil * phi + [15 / 12; -1 / 12; zeros(n - 4, 1);...
    -1 / 12; 15 / 12]) + 2 * phi .* (1 - phi) .* (2 * phi - 1 - (3 * L / (H * Tm)) * (T - Tm))) / tau;
phi = zeros(n, m);
%upper scaleth is lava, so height/scale=height0
phi(1 : int64(n - (n / scaling)), :) = ones(int64(n - (n / scaling)), m);
T0 = 1200+273.15;
T = T0 * ones(n, m);
T(1 : int64(n - (n / scaling)), :) = groundtemp * ones(int64(n - (n / scaling)), m);
% T = groundtemp : (T0 - 20) / (n - 1) : T0;
% T = T.' .* ones(1, m); 

%phase = zeros(length(xq) * int64(startflow / dt), 2);
const = alpha * cp;
% A = (height0 + height1) ^ 2;
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
while tq(k) < t1
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
                title(ax1, ['Days elapsed: ' num2str(t / (60 * 60 * 24))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                % Bottom plot
                ax2 = nexttile;
                plot(ax2, Ti)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [K]') 
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
                title(ax1, ['Days elapsed: ' num2str(t / (60 * 60 * 24))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                ax2 = nexttile;
                plot(ax2, Ti)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [K]') 
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
                title(ax1, ['Days elapsed: ' num2str(t / (60 * 60 * 24))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                % Bottom plot
                ax2 = nexttile;
                plot(ax2, Ti)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [K]') 
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
                title(ax1, ['Days elapsed: ' num2str(t / (60 * 60 * 24))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                ax2 = nexttile;
                plot(ax2, Ti)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [K]') 
                ylim([0 T0 * 1.1])             
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))     
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

%implement second layer
disp('Emplacing second layer.')
dn = int64(height1 / dx);
n = n + dn;
D = sparse(1 : n, 1 : n, -30 * ones(1, n) / 12, n, n);
E = sparse(2 : n, 1 : n - 1, 16 * ones(1, n - 1) / 12, n, n);
F = sparse(3 : n, 1 : n - 2, -ones(1, n - 2) / 12, n, n);
fdstencil = F + E + D + E' + F'; 

fphi = @(phi, T) (((omega / dx) ^ 2) * (fdstencil * phi + [15 / 12; -1 / 12; zeros(n - 4, 1);...
    -1 / 12; 15 / 12]) + 2 * phi .* (1 - phi) .* (2 * phi - 1 - (3 * L / (H * Tm)) * (T - Tm))) / tau;
phi = [phi; zeros(dn, m)];
T = [T; T0 * ones(dn, m)];
phii = phi(:, k);
Ti = T(:, k);
fphiq = zeros(n, 4);
fTq = zeros(n, 4);
dt = dt0;
fT = @(phi, Tl) (alpha / (dx ^ 2)) * (fdstencil * Tl + [15 * groundtemp / 12; -groundtemp / 12;...
    zeros(n - 4, 1); -TOne(Tl) / 12 ; ((16 * TOne(Tl)) - TTwo(Tl)) / 12]) + ((6 * L / cp) * phi .* (1 - phi) .* fphi(phi, Tl));
xq = dx * (0 : round((height / dx) / scale) : (height / dx));
i = 1;
startcounter = 1;
firsttime = 1;
while sum(Ti > Tend) > 0
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
            elseif qT > 2
                dt = min(2 * dt, 0.5 * (dx ^ 2) / alpha);
                Ti = Tf;
                startcounter = 1;
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
                title(ax1, ['Days elapsed: ' num2str(t / (60 * 60 * 24))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                % Bottom plot
                ax2 = nexttile;
                plot(ax2, Ti)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [K]') 
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
                title(ax1, ['Days elapsed: ' num2str(t / (60 * 60 * 24))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                ax2 = nexttile;
                plot(ax2, Ti)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [K]') 
                ylim([0 T0 * 1.1])             
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))        
            end
            k = k + 1;
        end
    else
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
            elseif (qphi > 2) && (qT > 2)
                dt = min(2 * dt, 0.5 * (dx ^ 2) / alpha);
                phii = phif;
                Ti = Tf;
                startcounter = 1;
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
                title(ax1, ['Days elapsed: ' num2str(t / (60 * 60 * 24))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                % Bottom plot
                ax2 = nexttile;
                plot(ax2, Ti)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [K]') 
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
                title(ax1, ['Days elapsed: ' num2str(t / (60 * 60 * 24))])
                xlabel(ax1,'height [m]')
                ylabel(ax1,'Φ (0 for liquid, 1 for solid)')
                ylim([-0.1 1 + 0.1]) 
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))            

                ax2 = nexttile;
                plot(ax2, Ti)
                xlabel(ax2,'height [m]')            
                ylabel(ax2,'T [K]') 
                ylim([0 T0 * 1.1])             
                xticks(0 : n / scale : n)
                xticklabels(num2str(xq'))        
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

elapsedTime = toc;
disp(['The MOL part of the simulation completed in ',...
    num2str(toc / 60), ' minutes.'])

% parse through data separately; back-track to find crossings
Tselect = T(int64((height - height0) / dx) + 1 : int64(height / dx), :);
phiselect = phi(int64((height - height0) / dx) + 1 : int64(height / dx), :);
xvec = dx * (1 : ((height + height1) / dx));
xsmallvec = dx * (1 : (height0 / dx));
flag1 = 0;
flag2 = 0;
flag3 = 0;
for i = flip(1 : k)
    if ~((sum(Tselect(:, i) > Tend) == 0) || flag1)
        tneeded = mean([tq(i - 1) tq(i)]);
        [~, I] = max(Tselect(:, i));
        xneeded = xsmallvec(min(length(xsmallvec), I));
        dataneeded = [tneeded xneeded];
        flag1 = 1;
    end
    % check for which point was the last to solidify
    if ~((sum(T(:, i) > Tm) == 0) || flag2)
        tlastoverall = mean([tq(i - 1) tq(i)]);
        [~, I] = min(phi(:, i));
        xlastoverall = xvec(I) - xvec(int64((height - height0) / dx));
        flag2 = 1;
        datalastoverall = [tlastoverall xlastoverall];
    end
    if ~((sum(Tselect(:, i) > Tm) == 0) || flag3)
        tlastfirst = mean([tq(i - 1) tq(i)]);
        [~, I] = min(phiselect(:, i));
        xlastfirst = xsmallvec(min(length(xsmallvec), I));
        flag3 = 1;
        datalastfirst = [tlastfirst xlastfirst];
    end
    if flag1 && flag2 && flag3 
        break
    end
end

heights = [height height0 height1];
times = [t1 t];
phi = phi(:, 1 : k);
T = T(:, 1 : k);
all = {d, phi, T, dx, k, heights, times, elapsedTime, dataneeded, datalastfirst, datalastoverall};