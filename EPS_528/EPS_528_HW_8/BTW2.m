Lq = 2 .^ (3 : 6);
nmax = 1e4;

data1 = BTW1(Lq(1), nmax, 0);
stuff = data1.occurrence;
count1 = stuff(1, :);
P1 = stuff(2, :) / length(data1.sizes);

data2 = BTW1(Lq(2), nmax, 0);
stuff = data2.occurrence;
count2 = stuff(1, :);
P2 = stuff(2, :) / length(data2.sizes);

data3 = BTW1(Lq(3), nmax, 0);
stuff = data3.occurrence;
count3 = stuff(1, :);
P3 = stuff(2, :) / length(data3.sizes);

data4 = BTW1(Lq(4), nmax, 0);
stuff = data4.occurrence;
count4 = stuff(1, :);
P4 = stuff(2, :) / length(data4.sizes);

figure(1)
loglog(count1, P1, count2, P2, count3, P3, count4, P4)