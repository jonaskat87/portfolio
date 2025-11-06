close all 

movienumber = 7; % change for movie name

% physical parameters
T = 8 * 60 * 60; % [s]
theta = 60; % latitude, in degrees
phi = 90; % local angle about vertical in degrees
% x goes along the channel, y goes across the channel
a = 6360e3; % [m]
g = 9.80665; % [m/s^2]
Omega = 7.2921e-5; % [rad/s]
H = 50; % [m]
D = 200e3; % [m]
L = 1000e3; % [m]
s = 5; % scaling parameter (>2)

% numerical parameters
dt = 2;
dx = L * 5e-3;
dy = D * 5e-3;
fluxlimiter = 'superbee'; % set flux limiter
% options for flux limiter include minmod, superbee, ospre,
%  vanleer, vanalbada

% simulation parameters
M = round(T / dt);
Nx = 2 * round(L / (2 * dx)) + 1;
Ny = 2 * round(D / (2 * dy)) + 1; 
t = 0;

% initial conditions
soli = zeros(Nx, Ny, 3);
% % dam-breaking
% soli(1 : (Nx - 1) / 2, :, 1) = H + (H / s);
% soli((Nx + 1) / 2 : end, :, 1) = H - (H / s);
% water slam
soli(:, :, 1) = H;
soli(1 : (Nx - 1) / 2, :, 2) = H;
soli((Nx + 1) / 2 : end, :, 2) = -H;
% % hump of water
% xq = 1 : Nx;
% soli(:, :, 1) = (H + (H / s) * (exp(-(400 *...
% ((xq / Nx) - 1 / 2).^2)) - 1)).' * ones(1, Ny);

% plotting parameters
plotno = Inf;
plotevery = 20;
plotno = min(min(Nx, plotno), Ny);
xvec = 1 : round(Nx / plotno) : Nx;
yvec = 1 : round(Ny / plotno) : Ny;
[X, Y] = meshgrid(dx * xvec - L / 2, dy * yvec - D / 2);

f0 = 2 * Omega * sin(theta * pi / 180);
beta = 2 * Omega * cos(theta * pi / 180) / a;
xq = 1 : Nx;
yq = 1 : Ny;
[Xq, Yq] = meshgrid(dx * xq - L / 2, dy * yq - D / 2);
Xq = Xq.';
Yq = Yq.';
f = f0 + beta * (Xq * sin(phi * pi / 180) + Yq * cos(phi * pi / 180));

moviename = ['MUSCL_', num2str(movienumber), '.mp4'];

fh = figure(1);
fh.WindowState = 'maximized';
vidfile = VideoWriter(moviename, 'MPEG-4');
vidfile.FrameRate = 5;
open(vidfile);
for k = 1 : M    
    solf = func(soli, f, dt, dx, dy, Nx, Ny, g, fluxlimiter);
    t = t + dt;
    
    % plotting
    if round(k / plotevery) == (k / plotevery)
        % height
        tiledlayout(2, 3) 
        ax1 = nexttile;
        hold on
            h = surf(ax1, X / 1e3, Y / 1e3, solf(xvec, yvec, 1).',...
                'FaceColor', 'interp', 'FaceAlpha',0.8);
            ground = surf(ax1, X / 1e3, Y / 1e3,...
                zeros(length(yvec), length(xvec)),...
                'FaceColor', [0.8 0.8 0.2]);
        hold off
        set(h, 'LineStyle', 'none');
        set(ground, 'LineStyle', 'none');
        title(ax1, 'Water height', 'FontWeight', 'normal')
        xlabel('x [km]')
        ylabel('y [km]')
        zlabel('h [m]')
        zlim([0 - (H / s) H + 2 * (H / s)])
        grid on
        colorbar('southoutside')
        caxis([H / 2 3 * H / 2])
        view(ax1, -phi, 30)
        
        % x-velocity 
        ax2 = nexttile;
        h = surf(ax2, X / 1e3, Y / 1e3,...
            (solf(xvec, yvec, 2) ./ solf(xvec, yvec, 1)).',...
            'FaceColor', 'interp','FaceAlpha',0.8);
        set(h, 'LineStyle', 'none');
        [title1, title2] = title(ax2,...
            ['Minutes elapsed: ' num2str(t / 60)],...
            'Velocity along the channel');
        title1.FontSize = 16;
        xlabel('x [km]')
        ylabel('y [km]')
        zlabel('u (x-velocity) [m/s]')
        zlim([-1 1] * 40 / (s - 1))
        grid on
        colorbar('southoutside')
        caxis([-1 1] * 40 / (s - 1))
        view(ax2, -phi, 30)
        
        % y-velocity
        ax3 = nexttile;
        h = surf(ax3, X / 1e3, Y / 1e3,...
            (solf(xvec, yvec, 3) ./ solf(xvec, yvec, 1)).',...
            'FaceColor', 'interp','FaceAlpha',0.8);
        set(h, 'LineStyle', 'none');
        title(ax3, 'Velocity across the channel',...
            'FontWeight', 'normal')
        xlabel('x [km]')
        ylabel('y [km]')
        zlabel('v (y-velocity) [m/s]')
        zlim([-1 1] * 15 / (s - 1))
        grid on
        colorbar('southoutside')
        caxis([-1 1] * 15 / (s - 1))
        view(ax3, -phi, 30)
        
        ax4 = nexttile;
        hold on
            h = surf(ax4, X / 1e3, Y / 1e3, solf(xvec, yvec, 1).',...
                'FaceColor', 'interp','FaceAlpha',0.8);
            ground = surf(ax4, X / 1e3, Y / 1e3,...
                zeros(length(yvec), length(xvec)),...
                'FaceColor', [0.8 0.8 0.2]);
        hold off
        set(h, 'LineStyle', 'none');
        set(ground, 'LineStyle', 'none');
        xlabel('x [km]')
        ylabel('y [km]')
        zlabel('h (water height) [m]')
        zlim([0 H + 2 * (H / s)])
        grid on
        caxis([H / 2 3 * H / 2])
        view(ax4, 90 - phi, 30)
        
        ax5 = nexttile;
        h = surf(ax5, X / 1e3, Y / 1e3,...
            (solf(xvec, yvec, 2) ./ solf(xvec, yvec, 1)).',...
            'FaceColor', 'interp','FaceAlpha',0.8);
        set(h, 'LineStyle', 'none');
        xlabel('x [km]')
        ylabel('y [km]')
        zlabel('u (x-velocity) [m/s]')
        zlim([-1 1] * 40 / (s - 1))
        grid on
        caxis([-1 1] * 40 / (s - 1))
        view(ax5, 90 - phi, 30)
        
        ax6 = nexttile;
        h = surf(ax6, X / 1e3, Y / 1e3,...
            (solf(xvec, yvec, 3) ./ solf(xvec, yvec, 1)).',...
            'FaceColor', 'interp','FaceAlpha',0.8);
        set(h, 'LineStyle', 'none');
        xlabel('x [km]')
        ylabel('y [km]')
        zlabel('v (y-velocity) [m/s]')
        zlim([-1 1] * 15 / (s - 1))
        grid on
        caxis([-1 1] * 15 / (s - 1))
        view(ax6, 90 - phi, 30)
        
        colormap(ax1, cold)
        colormap(ax2, jet)
        colormap(ax3, jet)
        colormap(ax4, cold)
        colormap(ax5, jet)
        colormap(ax6, jet)
        
    	drawnow
        F = getframe(fh); 
        writeVideo(vidfile, F);
    end   
    soli = solf;
end
close(vidfile)

function F = xflux(sol, Nx, Ny, grav)
    F = zeros(Nx, Ny, 3);
    ETA = sol(:, :, 1);
    U = sol(:, :, 2) ./ ETA;
    V = sol(:, :, 3) ./ ETA;
    F(:, :, 1) = U;
    F(:, :, 2) = U .* U + 0.5 * grav * ETA;
    F(:, :, 3) = U .* V;
    F = ETA .* F;
end

function G = yflux(sol, Nx, Ny, grav)
    G = zeros(Nx, Ny, 3);
    ETA = sol(:, :, 1);
    U = sol(:, :, 2) ./ ETA;
    V = sol(:, :, 3) ./ ETA;
    G(:, :, 1) = V;
    G(:, :, 2) = U .* V;
    G(:, :, 3) = V .* V + 0.5 * grav * ETA;
    G = ETA .* G;
end

function Fspeed = aF(solL, solR, g)
    ETAL = solL(:, :, 1);
    UL = solL(:, :, 2) ./ ETAL;
    ETAR = solR(:, :, 1);
    UR = solR(:, :, 2) ./ ETAR;
    FcompareL = repmat(UL, [1 1 3]);
    FcompareL(:, :, 2) = FcompareL(:, :, 2) - sqrt(g * ETAL);
    FcompareL(:, :, 3) = FcompareL(:, :, 3) + sqrt(g * ETAL);
    FcompareR = repmat(UR, [1 1 3]);
    FcompareR(:, :, 2) = FcompareR(:, :, 2) - sqrt(g * ETAR);
    FcompareR(:, :, 3) = FcompareR(:, :, 3) + sqrt(g * ETAR);
    largestL = max(abs(FcompareL), [], 3);
    largestR = max(abs(FcompareR), [], 3);
    Fspeed = max(largestL, largestR);
end

function Gspeed = aG(solL, solR, g)
    ETAL = solL(:, :, 1);
    VL = solL(:, :, 3) ./ ETAL;
    ETAR = solR(:, :, 1);
    VR = solL(:, :, 3) ./ ETAR;
    GcompareL = repmat(VL, [1 1 3]);
    GcompareL(:, :, 2) = GcompareL(:, :, 2) - sqrt(g * ETAL);
    GcompareL(:, :, 3) = GcompareL(:, :, 3) + sqrt(g * ETAL);
    GcompareR = repmat(VR, [1 1 3]);
    GcompareR(:, :, 2) = GcompareR(:, :, 2) - sqrt(g * ETAR);
    GcompareR(:, :, 3) = GcompareR(:, :, 3) + sqrt(g * ETAR);
    largestL = max(abs(GcompareL), [], 3);
    largestR = max(abs(GcompareR), [], 3);
    Gspeed = max(largestL, largestR);
end

function phi = fluxlim(vec, fluxlimiter)
    if strcmp(fluxlimiter, 'minmod')
    	phi = max(0, min(1, vec)); % minmod
    elseif strcmp(fluxlimiter, 'ospre')
        phi = 1.5 * ((vec .* vec) + vec) ./...
            ((vec .* vec) + vec + 1); % ospre
        phi(isnan(phi)) = 1.5; 
    elseif strcmp(fluxlimiter, 'superbee')
        phi = zeros(length(vec(:, 1, 1)),...
            length(vec(1, :, 1)), 3, 3); % superbee
        phi(:, :, :, 2) = min(2 * vec, 1);
        phi(:, :, :, 3) = min(vec, 2);
        phi = max(phi, [], 4);
    elseif strcmp(fluxlimiter, 'vanleer')
        phi = (vec + abs(vec)) ./ (1 + abs(vec)); % vanleer
        phi(isnan(phi)) = 2;       
    elseif strcmp(fluxlimiter, 'vanalbada')
        phi = ((vec .* vec) + vec) ./...
            ((vec .* vec ) + 1); % vanalbada(1)
        phi(isnan(phi)) = 1;    
    else
        phi = zeros(length(vec(:, 1, 1)),...
            length(vec(1, :, 1)), 3, 3); % superbee (default)
        phi(:, :, :, 2) = min(2 * vec, 1);
        phi(:, :, :, 3) = min(vec, 2);
        phi = max(phi, [], 4);    
    end
end

function dF = dFfluxdiff(soli, Nx, Ny, g, fluxlimiter)
    % transform into u and v for BCs
    soli(:, :, 2 : 3) = soli(:, :, 2 : 3) ./ soli(:, :, 1);
    
    % implement Neumann (no-penetration) BCs
    % x = -D / 2
    solba = circshift(soli, 1, 1);
    solba(1, :, :) = (4 * soli(1, :, :) - soli(2, :, :)) / 3;
    
    % x = -D / 2 - 
    solbaba = circshift(solba, 1, 1);
    solbaba(1, :, :) = soli(1, :, :); 

    % x = D / 2
    solfr = circshift(soli, -1, 1);
    solfr(end, :, :) = (4 * soli(end, :, :)...
        - soli(end - 1, :, :)) / 3;

    % x = D / 2 +
    solfrfr = circshift(solfr, -1, 1);
    solfrfr(end, :, :) = soli(end, :, :);
    
    % transform back into uh and vh
    soli(:, :, 2 : 3) = soli(:, :, 2 : 3) .* soli(:, :, 1);
    solba(:, :, 2 : 3) = solba(:, :, 2 : 3) .* solba(:, :, 1);
    solbaba(:, :, 2 : 3) = solbaba(:, :, 2 : 3) .* solbaba(:, :, 1);
    solfr(:, :, 2 : 3) = solfr(:, :, 2 : 3) .* solfr(:, :, 1);
    solfrfr(:, :, 2 : 3) = solfrfr(:, :, 2 : 3) .* solfrfr(:, :, 1);
    
    soliba = soli - solba;
    solfri = solfr - soli;
    solfrfrfr = solfrfr - solfr;
    oneflux = fluxlim(soliba ./ solfri, fluxlimiter);
    solLplus = soli + 0.5 * oneflux .* solfri;
    solLminus = solba + 0.5 * fluxlim((solba - solbaba) ./...
        soliba, fluxlimiter) .* soliba; 
    solRplus = solfr - 0.5 * fluxlim(solfri ./...
        solfrfrfr, fluxlimiter) .* solfrfrfr;
    solRminus = soli - 0.5 * oneflux .* solfri; 
    
    Fminus = 0.5 * (xflux(solRminus, Nx, Ny, g) +...
        xflux(solLminus, Nx, Ny, g) - aF(solLminus, solRminus, g)...
        .* (solRminus - solLminus));
    Fplus = 0.5 * (xflux(solRplus, Nx, Ny, g) +...
        xflux(solLplus, Nx, Ny, g) - aF(solLplus, solRplus, g)...
        .* (solRplus - solLplus));
    dF = Fplus - Fminus;
end

function dG = dGfluxdiff(soli, Nx, Ny, g, fluxlimiter)
    % transform into u and v for BCs
    soli(:, :, 2 : 3) = soli(:, :, 2 : 3) ./ soli(:, :, 1);

    % implement Neumann (no-penetration) BCs, aside for v at the walls,
    % which will have kinematic (homogeneous Dirichlet) BCs, i.e., v=0
    % y = 0 
    solbot = circshift(soli, 1, 2); 
    solbot(:, 1, 1 : 2) = (4 * soli(:, 1, 1 : 2)...
        - soli(:, 2, 1 : 2)) / 3;
    solbot(:, 1, 3) = 0;
    
    % y = 0 - 
    solbotbot = circshift(solbot, 1, 2);
    solbotbot(:, 1, 1 : 2) = soli(:, 1, 1 : 2);
    solbotbot(:, 1, 3) = 0; 
    
    % y = L 
    soltop = circshift(soli, -1, 2);
    soltop(:, end, 1 : 2) = ...
        (4 * soli(:, end, 1 : 2) - soli(:, end - 1, 1 : 2)) / 3;
    soltop(:, end, 3) = 0;
    
    % y = L +
    soltoptop = circshift(soltop, -1, 2);
    soltoptop(:, end, 1 : 2) = soli(:, end, 1 : 2);
    soltoptop(:, end, 3) = 0;
    
    % transform back into uh and vh
    soli(:, :, 2 : 3) = soli(:, :, 2 : 3) .* soli(:, :, 1);
    solbot(:, :, 2 : 3) = solbot(:, :, 2 : 3) .* solbot(:, :, 1);
    solbotbot(:, :, 2 : 3) = solbotbot(:, :, 2 : 3)...
        .* solbotbot(:, :, 1);
    soltop(:, :, 2 : 3) = soltop(:, :, 2 : 3) .* soltop(:, :, 1);
    soltoptop(:, :, 2 : 3) = soltoptop(:, :, 2 : 3)...
        .* soltoptop(:, :, 1);
    
    solibot = soli - solbot;
    soltopi = soltop - soli;
    soltoptoptop = soltoptop - soltop;
    oneflux = fluxlim(solibot ./...
        (soltop - soli), fluxlimiter);
    solLplus = soli + 0.5 * oneflux .* soltopi;
    solLminus = solbot + 0.5 * fluxlim((solbot - solbotbot) ./...
        solibot, fluxlimiter) .* solibot; 
    solRplus = soltop - 0.5 * fluxlim(soltopi ./...
        soltoptoptop, fluxlimiter) .* soltoptoptop;
    solRminus = soli - 0.5 * oneflux .* soltopi; 
    
    Gminus = 0.5 * (yflux(solRminus, Nx, Ny, g) +...
        yflux(solLminus, Nx, Ny, g) - aG(solLminus, solRminus, g)...
        .* (solRminus - solLminus));
    Gplus = 0.5 * (yflux(solRplus, Nx, Ny, g) +...
        yflux(solLplus, Nx, Ny, g) - aG(solLplus, solRplus, g)...
        .* (solRplus - solLplus));
    dG = Gplus - Gminus;
end

function vec = cor(sol, f, Nx, Ny)
    vec = zeros(Nx, Ny, 2);
    vec(:, :, 1) = f .* sol(:, :, 3);
    vec(:, :, 2) = -f .* sol(:, :, 2);
end

function solf = func(sol, f, dt, dx, dy, Nx, Ny, g, fluxlimiter)
% Strang splitting with midpoint method for all timesteps
    dt = dt / 2;
    k1 = -((dFfluxdiff(sol, Nx, Ny, g, fluxlimiter) / dx) + ...
        (dGfluxdiff(sol, Nx, Ny, g, fluxlimiter) / dy));
    k2 = -((dFfluxdiff(sol + dt * k1 / 2,...
        Nx, Ny, g, fluxlimiter) / dx) + ...
        (dGfluxdiff(sol + dt * k1 / 2, Nx, Ny, g, fluxlimiter) / dy));
    solf = sol + dt * k2; % first time step
    
    dt = 2 * dt;
    % do not need to update h for Coriolis terms
    k1(:, :, 2 : 3) = cor(solf, f, Nx, Ny);
    k2 = cor(solf + dt * k1 / 2, f, Nx, Ny);
    solf(:, :, 2 : 3) = solf(:, :, 2 : 3) + dt * k2; % second time step
    
    dt = dt / 2;
    k1 = -((dFfluxdiff(solf, Nx, Ny, g, fluxlimiter) / dx) + ...
        (dGfluxdiff(solf, Nx, Ny, g, fluxlimiter) / dy));
    k2 = -((dFfluxdiff(solf + dt * k1 / 2,...
        Nx, Ny, g, fluxlimiter) / dx) + ...
        (dGfluxdiff(solf + dt * k1 / 2, Nx, Ny, g, fluxlimiter) / dy));
    solf = solf + dt * k2; % final time step
end
