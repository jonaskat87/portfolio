% initial conditions 
tol = 1e-10;
t0 = 0; % initial time
tf = 50; % final time
tspan = [t0 tf];

Iq = 43 : 51; 
V0 = [-50 50 -40 50 50 -30 -50]; 
n0 = [0.1 0.1 0.6 0.5 0.9 0.8 0.35];
Vq = -75 : 0.5 : 50; % domain values for nullclines

C = 1;
gL = 1;
EL = -78;
gNa = 4;
ENa = 60;
Vm = -30;
km = 7;
gK = 4;
EK = -90;
Vn = -45;
kn = 5;
tauV = 1;

for I = Iq
    simulation(tol, V0, n0, tspan, Vq,...
    I, C, gL, EL, gNa, ENa, Vm, km, gK, EK, Vn, kn, tauV)
end