tic
Lq = 2 .^ (3 : 6);
nmax = 1e5;

data1 = BTWfirst(Lq(1), nmax, 0, 0);
stuff = data1.occurrence;
P1 = stuff(1, :) / nmax;
domain1 = stuff(2, :);
Pshorter = P1(2 : ceil(end / 5));
domainshorter = domain1(2 : ceil(end / 5));
Pnonzero = Pshorter(intersect(...
    intersect(find(domainshorter < 21),...
    find(Pshorter)), find(domainshorter)));
domainnonzero = domainshorter(intersect(intersect(find(...
    domainshorter < 21), find(...
    Pshorter)), find(domainshorter)));
m1 = sum((log(Pnonzero) - mean(log(Pnonzero)))...
    .* (log(domainnonzero) - mean(log(domainnonzero))))...
    / sum((log(domainnonzero) -...
    mean(log(domainnonzero))) .^ 2);

data2 = BTWfirst(Lq(2), nmax, 0, 0);
stuff = data2.occurrence;
P2 = stuff(1, :) / nmax;
domain2 = stuff(2, :);
Pshorter = P2(2 : ceil(end / 5));
domainshorter = domain2(2 : ceil(end / 5));
Pnonzero = Pshorter(intersect(...
    intersect(find(domainshorter < 21),...
    find(Pshorter)), find(domainshorter)));
domainnonzero = domainshorter(intersect(intersect(find(...
    domainshorter < 21), find(...
    Pshorter)), find(domainshorter)));
m2 = sum((log(Pnonzero) - mean(log(Pnonzero)))...
    .* (log(domainnonzero) - mean(log(domainnonzero))))...
    / sum((log(domainnonzero) -...
    mean(log(domainnonzero))) .^ 2);

data3 = BTWfirst(Lq(3), nmax, 0, 0);
stuff = data3.occurrence;
P3 = stuff(1, :) / nmax;
domain3 = stuff(2, :);
Pshorter = P3(2 : ceil(end / 5));
domainshorter = domain3(2 : ceil(end / 5));
Pnonzero = Pshorter(intersect(...
    intersect(find(domainshorter < 21),...
    find(Pshorter)), find(domainshorter)));
domainnonzero = domainshorter(intersect(intersect(find(...
    domainshorter < 21), find(...
    Pshorter)), find(domainshorter)));
m3 = sum((log(Pnonzero) - mean(log(Pnonzero)))...
    .* (log(domainnonzero) - mean(log(domainnonzero))))...
    / sum((log(domainnonzero) -...
    mean(log(domainnonzero))) .^ 2);

data4 = BTWfirst(Lq(4), nmax, 0, 0);
stuff = data4.occurrence;
P4 = stuff(1, :) / nmax;
domain4 = stuff(2, :);
Pshorter = P4(2 : ceil(end / 5));
domainshorter = domain4(2 : ceil(end / 5));
Pnonzero = Pshorter(intersect(...
    intersect(find(domainshorter < 21),...
    find(Pshorter)), find(domainshorter)));
domainnonzero = domainshorter(intersect(intersect(find(...
    domainshorter < 21), find(...
    Pshorter)), find(domainshorter)));
m4 = sum((log(Pnonzero) - mean(log(Pnonzero)))...
    .* (log(domainnonzero) - mean(log(domainnonzero))))...
    / sum((log(domainnonzero) -...
    mean(log(domainnonzero))) .^ 2);

figure(1)
loglog(domain1, P1, domain2, P2, domain3, P3, domain4, P4)
title('Avalanche probability vs. size', 'interpreter',...
    'latex', 'fontsize', 18)
xlabel('Avalanche size ($s$)', 'interpreter',...
    'latex', 'fontsize', 14)
ylabel('Avalanche probability ($P(s; L)$)',...
    'interpreter', 'latex', 'fontsize', 14)
legend(['L=8 (slope: ', num2str(m1), ')'],...
    ['L=16 (slope: ', num2str(m2), ')'],...
    ['L=32 (slope: ', num2str(m3), ')'],...
    ['L=64 (slope: ', num2str(m4), ')'], 'interpreter', 'latex')
toc