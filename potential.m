clear all
close all

n = 10e-9; %concentration
a = 1000; %nm
L = 200; % deplent size nm
phi = n*4.0/3.0*pi*L^3;
att = 3.5;
rm=2.4;
r = 0:0.01:rm;
r = r*a;
kappa = 50;
Bpp = 2.29*a;
u_rep = Bpp*exp(-kappa*(r -2*a)/a);

P=n*(1+phi+phi^2-phi^3)/(1-phi)^3;
P = 5.8e-8; % P has unit of KT nm^3
u_ao = -P*pi*(4.0/3.0*(a+L)^3.*(1-3*r./(4*(a+L))+r.^3./(16*(a+L)^3)));

figure(1)
plot(r/a,u_rep);
xlabel('r/a');
ylabel('u/kT');
xlim([2, 2+0.4])
ylim([-25 10]);
figure(2)
plot(r,u_ao);
figure(1)
hold on
plot(r/a,u_rep+u_ao);
xlabel('r/a');
ylabel('u/kT');
xlim([2, 2+0.4])
ylim([-25 5]);
P = 8.8e-8; % P has unit of KT nm^3
u_ao = -P*pi*(4.0/3.0*(a+L)^3.*(1-3*r./(4*(a+L))+r.^3./(16*(a+L)^3)));
plot(r/a,u_rep+u_ao);
P = 12.8e-8; % P has unit of KT nm^3
u_ao = -P*pi*(4.0/3.0*(a+L)^3.*(1-3*r./(4*(a+L))+r.^3./(16*(a+L)^3)));
plot(r/a,u_rep+u_ao);
P = 16.8e-8; % P has unit of KT nm^3
u_ao = -P*pi*(4.0/3.0*(a+L)^3.*(1-3*r./(4*(a+L))+r.^3./(16*(a+L)^3)));
plot(r/a,u_rep+u_ao);
P = 20.8e-8; % P has unit of KT nm^3
u_ao = -P*pi*(4.0/3.0*(a+L)^3.*(1-3*r./(4*(a+L))+r.^3./(16*(a+L)^3)));
plot(r/a,u_rep+u_ao);

legend('DL repulsion','DL+AO (P=5.8e-8)','DL+AO (P=8.8e-8)','DL+AO (P=12.8e-8)','DL+AO (P=16.8e-8)','DL+AO (P=20.8e-8)')