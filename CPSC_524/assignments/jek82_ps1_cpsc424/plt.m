% plots output of triad (MFlops vs N; problem 2)
data = readmatrix('tmp.txt');
f = figure('visible','off');
semilogx(data(:,1), data(:,2));
title("Benchmarking the vector triad kernel");
xlabel("N (array length; log scale)");
ylabel("MFLOPS (linear scale)");
grid on
saveas(f,'plt','pdf');
