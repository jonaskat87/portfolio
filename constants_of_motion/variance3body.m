function [sigmas, choices] = variance3body(k, u0, ti, tf,...
    M, G, tol, snr, sketch, choose)
% k determines the number of constants of motion we want to look at
% M = [1 0.9 1];
% G = 1 ;
% u0 = [1 -0.3 1 1 2 0 0 0 0 0 0 0];
% to = 0;
% tf = 200;
% tol = 1e-8;
% snr = 50; % signal-to-noise ratio
% sketch = 0; % do you want to use a low-rank matrix sketch or not?
% choose = 1,2,3,4,5; % which constant of motion do you want to consider?

m1 = M(1) ;
m2 = M(2) ;
m3 = M(3) ;
f = @(t, Y) [Y(7);Y(8);Y(9);Y(10);Y(11);Y(12);...
    G*(m2*(Y(3)-Y(1))/((Y(3)-Y(1))^2+(Y(4)-Y(2))^2)^(3/2)...
    +m3*(Y(5)-Y(1))/((Y(5)-Y(1))^2+(Y(6)-Y(2))^2)^(3/2));...
    G*(m2*(Y(4)-Y(2))/((Y(3)-Y(1))^2+(Y(4)-Y(2))^2)^(3/2)...
    +m3*(Y(6)-Y(2))/((Y(5)-Y(1))^2+(Y(6)-Y(2))^2)^(3/2));...
    G*(m1*(Y(1)-Y(3))/((Y(3)-Y(1))^2+(Y(4)-Y(2))^2)^(3/2)...
    +m3*(Y(5)-Y(3))/((Y(5)-Y(3))^2+(Y(6)-Y(4))^2)^(3/2));...
    G*(m1*(Y(2)-Y(4))/((Y(3)-Y(1))^2+(Y(4)-Y(2))^2)^(3/2)...
    +m3*(Y(6)-Y(4))/((Y(5)-Y(3))^2+(Y(6)-Y(4))^2)^(3/2));...
    G*(m1*(Y(1)-Y(5))/((Y(1)-Y(5))^2+(Y(6)-Y(2))^2)^(3/2)...
    +m2*(Y(3)-Y(5))/((Y(5)-Y(3))^2+(Y(6)-Y(4))^2)^(3/2));...
    G*(m1*(Y(2)-Y(6))/((Y(5)-Y(1))^2+(Y(6)-Y(2))^2)^(3/2)...
    +m2*(Y(4)-Y(6))/((Y(5)-Y(3))^2+(Y(6)-Y(4))^2)^(3/2))];
options = odeset('RelTol', tol, 'AbsTol', tol);
[~, u] = ode45(f, [ti, tf], u0, options);
u = awgn(u, snr);

if choose == 1 % Hamiltonian
    g = @(U) [1 ./ sqrt((U(:, 3)-U(:, 1)).^2+(U(:, 4)-U(:, 2)).^2)...
    1 ./ sqrt((U(:, 5)-U(:, 1)).^2+(U(:, 6)-U(:, 2)).^2)...
    1 ./ sqrt((U(:, 5)-U(:, 3)).^2+(U(:, 6)-U(:, 4)).^2)...
    U(:, 7) .^ 2 U(:, 8) .^ 2 U(:, 9) .^ 2 ...
    U(:, 10) .^ 2 U(:, 11) .^ 2 U(:, 12) .^ 2];
elseif choose == 2 % x linear momentum
    g = @(U) [U(:, 7) U(:, 9) U(:, 11)];    
elseif choose == 3 % y linear momentum
    g = @(U) [U(:, 8) U(:, 10) U(:, 12)];
elseif choose == 4 % angular momentum 
    g = @(U) [U(:, 1) .* U(:, 8) - U(:, 2) .* U(:, 7)...
        U(:, 3) .* U(:, 10) - U(:, 4) .* U(:, 9)...
        U(:, 5) .* U(:, 12) - U(:, 6) .* U(:, 11)];    
elseif choose == 5 % all at once
    g = @(U) [1 ./ sqrt((U(:, 3)-U(:, 1)).^2+(U(:, 4)-U(:, 2)).^2)...
        1 ./ sqrt((U(:, 5)-U(:, 1)).^2+(U(:, 6)-U(:, 2)).^2)...
        1 ./ sqrt((U(:, 5)-U(:, 3)).^2+(U(:, 6)-U(:, 4)).^2)...
        U(:, 7) .^ 2 U(:, 8) .^ 2 U(:, 9) .^ 2 ...
        U(:, 10) .^ 2 U(:, 11) .^ 2 U(:, 12) .^ 2 ...
        U(:, 7 : 12) U(:, 1) .* U(:, 8) - U(:, 2) .* U(:, 7)...
        U(:, 3) .* U(:, 10) - U(:, 4) .* U(:, 9)...
        U(:, 5) .* U(:, 12) - U(:, 6) .* U(:, 11)];
else 
    disp('Illegal choice of constant of motion.')
    return
end
gmean = sum(g(u)) / length(u(:, 1));
gsize = g(ones(1, 12));

M = g(u);
M = M - gmean .* ones(size(M));
if sketch
%     [~, S, V] = svdsketch(M);
%     sigmas = diag(S);
%     sigmas = sigmas(end - k + 1 : end);
%     choices = V(:, end - k + 1 : end) ./...
%         max(abs(V(:, end - k + 1 : end)));
    [U, S, V] = svdsketch(M);
    sigmas = diag(S);
    Mapprox = zeros(size(M));
    for i = 1 : length(gsize) - k
        Mapprox = Mapprox + sigmas(i) * U(:, i) * transpose(V(:, i));
    end
    choices = null(Mapprox);
    sigmas = vecnorm(M * choices) / norm(M);
else
%     [U, S, V] = svds(M, length(gsize) - k);
%     sigmas = diag(S);
%     Mapprox = zeros(size(M));
%     for i = 1 : length(gsize) - k
%         Mapprox = Mapprox + sigmas(i) * U(:, i) * transpose(V(:, i));
%     end
%     choices = null(Mapprox);
%     sigmas = vecnorm(M * choices) / norm(M);
    [V, D] = eigs(transpose(M) * M, k, 'smallestabs');
    sigmas = sqrt(diag(D)) / norm(M);
    choices = V;
end