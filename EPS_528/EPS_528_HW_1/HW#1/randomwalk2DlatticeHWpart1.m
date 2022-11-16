function [ numberinfected ] = randomwalk2DlatticeHWpart1(N, M, Q, L, n)
%N = linear dimension (lattice size is N^2)
%M = total population number
%Q = number of initial infected people
%L = days to recover
%n = number of total iterations

%start with agents in random locations
loc = randi(N, 2, M);
infections = zeros(1, M); %to track infections
%zero is healthy, nonzero is infected
infections(randperm(M, Q)) = L; %set Q random agent(s) to be infected
numberinfected = zeros(1, n); %to track number infected
numberinfected(1) = sum(infections > 0);

%sample directions (to add on at each iteration)
dirq = [0 0 0 1 -1; 0 1 -1 0 0];
for i = 2 : n
    infections = max(0, infections - 1); %allow infection to pass 
    loc = mod(loc + dirq(:, randi(5, 1, M)), N);
    for j = find(infections > 0) % go through infected agents
        infloc = loc(:, j);
        % check where locations line up with infected
        equality = (loc == infloc); 
        % check which were infected (but remove already infected)
        infections((sum(equality) == 2) & (infections == 0)) = L; 
    end 
    numberinfected(i) = sum(infections > 0);
end
