function [ data ] = simulation2extraction(N, ind, Tend, actual)

cd 'C:\Users\jonas\OneDrive\Documents\MATLAB\newemplacementresults'
mq = [0.1 : 0.1 : 1 2 3 4 5 10 15 20];
height0 = mq(ind);
height1 = mq(ind);
h = mean([height0 height1]);
if actual
    fname = 'propconstant.mat';
    m = matfile(fname);
    all = m.(['case', num2str(ind)]);
    datalastoverall = all{11};
    th = datalastoverall(1);
else
    alpha = 3.75e-7;
    th = h ^ 2  / alpha;
end

fname = ['emplacementresults_', num2str(height0), '_', num2str(height1),...
 '_', num2str(round(Tend)), 'K.mat'];
m = matfile(fname);

x = zeros(1, N); %t_{emp}/t_{h}
% fig 1a: time to full solidification
fig1avec = zeros(1, N); %y-variable
% fig 1b: time to solidification of first lobe only
fig1bvec = zeros(1, N);
% fig 2: time to Tend for first lobe
fig2vec = zeros(1, N);
% fig 3a: depth where full solidification
fig3avec = zeros(1, N);
% fig 3b: depth where solidification of first lobe
fig3bvec = zeros(1, N);
% fig 4: depth where first lobe reached Tend
fig4vec = zeros(1, N);
for i = 1 : N
    all = m.(['all', num2str(i)]);
    times = all{7};
    dataneeded = all{9};
    datalastfirst = all{10};
    datalastoverall = all{11};
    x(i) = times(1) / th;
    fig1avec(i) = datalastoverall(1) / th;
    fig1bvec(i) = datalastfirst(1) / th;
    fig2vec(i) = dataneeded(1) / th;
    fig3avec(i) = datalastoverall(2) / h;
    fig3bvec(i) = datalastfirst(2) / h;
    fig4vec(i) = dataneeded(2) / h;
end

data = [x; fig1avec; fig1bvec; fig2vec; fig3avec; fig3bvec; fig4vec];

cd ..