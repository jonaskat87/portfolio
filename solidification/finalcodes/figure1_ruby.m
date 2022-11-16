clear,clc;
%close all;

n=8;
height0=10;
height1=10;
Tend=570+273;
factor = 1;
Tm = 1070; %degree C

load('lajolla.mat');
blue=[0 102 204]./255;

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
T=T-273; % convert temperature to degree C

x_range=height0*4/dx:1:height0*6/dx;
x_range_first=height0*4/dx:1:height0*5/dx; % first lobe

% identifying the times steps for the first-lobe only
first_index=find(tq<times(1));
if (n>2)
    first_index(end)=[];
    if (n>5)
        first_index(end)=[];
    end
end
phi_first=phi(x_range_first,:);

% calculate solid and melt thickness during the first-lobe only period
solid_frac_first=sum(phi_first(:,first_index),1)./length(x_range_first);
melt_thickness_first=height0*(1-solid_frac_first);

solid_frac=sum(phi(x_range,:),1)./length(x_range);
melt_thickness_total=2*height0*(1-solid_frac);

% compile the solid and melt thickness during the entire simulation
solid_frac(first_index)=solid_frac_first;
melt_thickness_total(first_index)=melt_thickness_first;


% redefine x-range for plotting temperature field
x_range=height0*2/dx:1:height0*6/dx;
T=T(x_range,:);
phi=phi(x_range,:);
xvec=dx*x_range;
time_vec=tq / (60 * 60 * 24 * 30);

if (n<3)
    trim_index=find(time_vec<=16);
elseif (n<7)
    trim_index=find(time_vec<=16);
    
else
    trim_index=find(time_vec<=51);
end

time_vec=time_vec(trim_index);
T=T(:,trim_index);

[XX, YY]=meshgrid(time_vec, xvec);

%% Top row
figure(1);clf;
subplot(211);
surf(time_vec,xvec,T);view([0,90]);shading interp; grid off;
axis tight; caxis([25 1200]);colormap((flipud(lajolla)));

hold on;
contour3((time_vec),xvec,T,'LevelList',1080,'Color','k','linewidth',2);
plot3([0,(time_vec(end))],[40, 40],[2000, 2000],'w-','linewidth',2);
ylabel('height [m]');
h=colorbar('southoutside','color','k');
set(h,'Xdir','reverse');

xlabel('time (months)');


%% bottom row
if (n==1)
    tq=[tq(1:2),times(1),tq(3:end)];
    time_vec=tq / (60 * 60 * 24 * 30);
    melt_thickness_total=[melt_thickness_total(1:2),melt_thickness_total(2)+10,melt_thickness_total(3:end)];
end

if (n==2)
    tq=[tq(1:3),times(1),tq(4:end)];
    time_vec=tq / (60 * 60 * 24 * 30);
    melt_thickness_total=[melt_thickness_total(1:3),melt_thickness_total(3)+10,melt_thickness_total(4:end)];
end


if (n<3)
    trim_index=find(time_vec<=16);
elseif (n<7)
    trim_index=find(time_vec<=16);
    
else
    trim_index=find(time_vec<=51);
end
time_vec=time_vec(trim_index);
melt_thickness_total=melt_thickness_total(trim_index);

figure(1);subplot(212);
plot(time_vec,melt_thickness_total,'-','color',blue,'linewidth',2);axis([0 51 0 20]);

xlabel('time (months)');
ylabel('Total melt thickness (m)');





%% extra code
% % find phi=0.5 index;
% [I,J]=find(phi<0.99 & phi >0.01);
% % find T=570 index;
% [It,Jt]=find(T<571 & T > 569);
%plot3(sqrt(time_vec(J)),xvec(I),2000*ones(length(I),1),'k.');
%plot3(sqrt(time_vec(Jt)),xvec(It),2000*ones(length(It),1),'b.');



% figure(2);
% subplot(211);
% surf(time_vec,xvec,T);view([0,90]);shading interp; grid off;
% axis tight; caxis([25 1300]);colormap((lajolla));
% hold on;
% contour3(time_vec,xvec,T,'LevelList',1080,'Color','k')
%
%
% solid_frac=sum(phi,1)./length(xvec);
% figure(2);subplot(212);
% plot(time_vec,solid_frac,'-');axis tight;





