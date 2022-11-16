function [ sol ] = problem3c1(x10, x20, t0, tf)
tspan = [t0 tf]; % set initial and end times
u0 = [x10 x20]; % set initial values for x_{1} and x_{2}
opts = odeset('RelTol', 1e-12); % set error tolerance
func = @(t, x) [x(2) ; x(1) * (1 - x(1) ^ 2)]; % function handle
sol = ode45(func, tspan, u0, opts); % use ode45 to solve ODE