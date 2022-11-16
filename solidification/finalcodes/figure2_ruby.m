clear,clc;
%close all;
figure(1);clf;

height0=0.1;
cc=hot(16);
for n=6


height1=height0;

Tend=570+273;
factor = 1;
Tm = 1070; %degree C

fname = ['../data/emplacementresults_', num2str(height0), '_', num2str(height1),...
    '_', num2str(round(Tend)), 'K.mat'];
m = matfile(fname, 'Writable', true);
all = m.(['all', num2str(n)]);
phi = all{2};
T = all{3};
dx = all{4};
k = all{5};
heights = all{6};
times = all{7};
t = times(2);
empt = (times(1) - (t / k))  / (60 * 60 * 24 *30);
tq = (t / k) * (0 : round(k / factor) - 1);
tq=tq./(60 * 60 * 24 * 30);
empt_time=times(1)./(60 * 60 * 24 * 30);

T=T-273; % convert temperature to degree C

x_range=height0*4/dx:1:height0*6/dx;
x_range=uint32(x_range);
T=T(x_range,:);
phi=phi(x_range,:);
xvec=dx*x_range;

solid_frac=sum(phi,1)./length(x_range);
figure(1);subplot(9,1,n);
plot(tq,solid_frac,'-','linewidth',2,'color',cc(n,:));
xlabel('time [months]');ylabel('solid fraction []');
titlename=sprintf('slab thickness = %d m,emplaced at t= %2.5f months',height0, empt_time);
title(titlename);


end
%figure(1);subplot(911);title('slab thickness = 10m');

for index=2:364
    figure(2);
    plot(phi(:,index),x_range,'-o');
    pause;
end