% initial conditions 
tol = 1e-10;
t0 = 0; % initial time
tf = 50; % final time
tspan = [t0 tf];

Iq = 3 : 0.5 : 6; 
V0 = [20 * ones(1, 3) -10 -30 -40 -60]; 
n0 = [0 0.5 0.9 0 0 1 1] ;
Vq = -75 : 0.5 : 20; % domain values for nullclines

C = 1;
gL = 8;
EL = -80;
gNa = 20;
ENa = 60;
Vm = -20;
km = 15;
gK = 10;
EK = -90;
Vn = -25;
kn = 5;
tauV = 0.14;

for I = Iq
    simulation(tol, V0, n0, tspan, Vq,...
    I, C, gL, EL, gNa, ENa, Vm, km, gK, EK, Vn, kn, tauV)
end