L = 128;
nmax = 1e4;
pq = 0.25 : 0.25 : 1;

data1 = BTWnonconsfirst(L, nmax, 0, 0, pq(1));
stuff = data1.occurrence;
P1 = stuff(1, :) / nmax;
domain1 = stuff(2, :);
disp('Completed p=0.25.')

data2 = BTWnonconsfirst(L, nmax, 0, 0, pq(2));
stuff = data2.occurrence;
P2 = stuff(1, :) / nmax;
domain2 = stuff(2, :);
disp('Completed p=0.5.')

data3 = BTWnonconsfirst(L, nmax, 0, 0, pq(3));
stuff = data3.occurrence;
P3 = stuff(1, :) / nmax;
domain3 = stuff(2, :);
disp('Completed p=0.75.')

data4 = BTWnonconsfirst(L, nmax, 0, 0, pq(4));
stuff = data4.occurrence;
P4 = stuff(1, :) / nmax;
domain4 = stuff(2, :);
disp('Completed p=1.')

% we already know that it will roughly lie in this range
Pshorter = P4(2 : ceil(end / 2));
domainshorter = domain4(2 : ceil(end / 2));
% only do lsq when avalanche size is under 30
Pnonzero = Pshorter(intersect(...
    intersect(find(domainshorter < 30),...
    find(Pshorter)), find(domainshorter)));
domainnonzero = domainshorter(intersect(intersect(find(...
    domainshorter < 30), find(Pshorter)), find(domainshorter)));

P = @(s, p) (factorial(4 * s) ./ (factorial(s)...
    .* factorial(3 * s))) .* (p .^ (s - 1)) .*...
    ((1 - p) .^ (3 * s + 1)) ./ (3 * s + 1);
P1theory = (dot(P1(2 : end), P(domain1(2 : end), pq(1) / 4))...
    / dot(P(domain1(2 : end), pq(1) / 4), P(domain1(2 : end),...
    pq(1) / 4))) * P(domain1, pq(1) / 4);
P2theory = (dot(P2(2 : end), P(domain2(2 : end), pq(2) / 4))...
    / dot(P(domain2(2 : end), pq(2) / 4), P(domain2(2 : end),...
    pq(2) / 4))) * P(domain2, pq(2) / 4);
P3theory = (dot(P3(2 : end), P(domain3(2 : end), pq(3) / 4))...
    / dot(P(domain3(2 : end), pq(3) / 4), P(domain3(2 : end),...
    pq(3) / 4))) * P(domain3, pq(3) / 4);
P4theory = (dot(Pnonzero, P(domainnonzero, pq(4) / 4))...
    / dot(P(domainnonzero, pq(4) / 4), P(domainnonzero,...
    pq(4) / 4))) * P(domain4, pq(4) / 4);

figure(1)
loglog(domain1, [P1; P1theory])
title('Avalanche probability vs. size ($p=0.25$)',...
    'interpreter', 'latex', 'fontsize', 18)
xlabel('Avalanche size ($s$)', 'interpreter',...
    'latex', 'fontsize', 14)
ylabel('Avalanche probability ($P(s; p)$)',...
    'interpreter', 'latex', 'fontsize', 14)
legend('Numerical', 'Theory', 'interpreter', 'latex')

figure(2)
loglog(domain2, [P2; P2theory])
title('Avalanche probability vs. size ($p=0.5$)',...
    'interpreter', 'latex', 'fontsize', 18)
xlabel('Avalanche size ($s$)', 'interpreter',...
    'latex', 'fontsize', 14)
ylabel('Avalanche probability ($P(s; p)$)',...
    'interpreter', 'latex', 'fontsize', 14)
legend('Numerical', 'Theory', 'interpreter', 'latex')

figure(3)
loglog(domain3, [P3; P3theory])
title('Avalanche probability vs. size ($p=0.75$)',...
    'interpreter', 'latex', 'fontsize', 18)
xlabel('Avalanche size ($s$)', 'interpreter',...
    'latex', 'fontsize', 14)
ylabel('Avalanche probability ($P(s; p)$)',...
    'interpreter', 'latex', 'fontsize', 14)
legend('Numerical', 'Theory', 'interpreter', 'latex')

figure(4)
loglog(domain4, [P4; P4theory])
title('Avalanche probability vs. size ($p=1$)',...
    'interpreter', 'latex', 'fontsize', 18)
xlabel('Avalanche size ($s$)', 'interpreter',...
    'latex', 'fontsize', 14)
ylabel('Avalanche probability ($P(s; p)$)',...
    'interpreter', 'latex', 'fontsize', 14)
legend('Numerical', 'Theory', 'interpreter', 'latex')

figure(5)
loglog(domain1, P1, domain2, P2,...
    domain3, P3, domain4, P4)
title('Avalanche probability vs. size (all cases)',...
    'interpreter', 'latex', 'fontsize', 18)
xlabel('Avalanche size ($s$)', 'interpreter',...
    'latex', 'fontsize', 14)
ylabel('Avalanche probability ($P(s; p)$)',...
    'interpreter', 'latex', 'fontsize', 14)
legend('p=0.25', 'p=0.5', 'p=0.75',...
    'p=1', 'interpreter', 'latex')