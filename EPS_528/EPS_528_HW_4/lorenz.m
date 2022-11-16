function [ sol ] = lorenz(sigma, r, b, t0, tf, x0, y0, z0)

opts = odeset('RelTol', 1e-12, 'AbsTol', 1e-10);
func = @(t, u) [sigma * (u(2) - u(1)); r * u(1) - u(2)...
    - u(1) * u(3); u(1) * u(2) - b * u(3)];
u0 = [x0; y0; z0];
tspan = [t0 tf];
sol = ode45(func, tspan, u0, opts);
end