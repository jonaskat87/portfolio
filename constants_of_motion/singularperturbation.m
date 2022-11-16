epsilon = 1e-10;
A = @(e) [1 + e 0; 0 1 - e];
E = @(e) [-e e; e e];
[~, S0, V0] = svd(A(epsilon));
[~, Sp, Vp] = svd(A(epsilon) + E(epsilon));
sigmas0 = diag(S0);
choices0 = V0;
sigmasp = diag(Sp);
choicesp = Vp;
disp(norm(A(epsilon) * V0, 'fro'))
disp(norm((A(epsilon) + E(epsilon)) * Vp, 'fro'))
