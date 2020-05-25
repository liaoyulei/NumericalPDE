%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve Sod's problem with init  %
% (\pho_l,v_l,p_l)=(1,0,1)       %
% (\pho_r,v_r,p_r)=(0.125,0,0.1) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x1 = -0.5; %init space
x2 = 0.5; %final space
CFL = 1; %dt^3/dx^5
T = 0.1; %final time
N = 400; %number of space interval
dx = (x2 - x1) ./ N;
dt = (CFL .* dx.^5).^(1/3);
xx = x1+dx./2: dx: x2;
[rho0, v0, p0] = ExSolu(xx, 0);
[rhoEXACT, vEXACT, pEXACT] = ExSolu(xx, T);
[rhoWENO, vWENO, pWENO] = NuSolu(rho0, v0, p0, dx, dt, T, 'WENO');
[rhoFWENO, vFWENO, pFWENO] = NuSolu(rho0, v0, p0, dx, dt, T, 'FWENO');
subplot(2, 2, 1);
plot(xx, rhoEXACT, xx, rhoWENO, xx, rhoFWENO);
axis([-0.5, 0.5, -0.1, 1.1]);
xlabel('x');
ylabel('\rho');
title('T=0.1, N=400, density');
subplot(2, 2, 2);
plot(xx, vEXACT, xx, vWENO, xx, vFWENO);
axis([-0.5, 0.5, -0.1, 1.1]);
xlabel('x');
ylabel('v');
title('T=0.1, N=400, velocity');
subplot(2, 2, 3);
plot(xx, pEXACT, xx, pWENO, xx, pFWENO);
axis([-0.5, 0.5, -0.1, 1.1]);
xlabel('x');
ylabel('p');
title('T=0.1, N=400, pressure');
subplot(2, 2, 4);
plot(0, 0, 0, 0, 0, 0);
legend('exact', 'WENO5', 'FWENO5', 'Location', 'southwest');
axis off;
print('-depsc', 'Sod1.eps');

subplot(2, 2, 1);
plot(xx, rhoEXACT, xx, rhoWENO, xx, rhoFWENO);
axis([-0.13, 0.01, 0.3, 1.1]);
xlabel('x');
ylabel('\rho');
title('Zoom-in view near the rarefaction wave');
subplot(2, 2, 2);
plot(xx, rhoEXACT, xx, rhoWENO, xx, rhoFWENO);
axis([0.06, 0.12, 0.2, 0.5]);
xlabel('x');
ylabel('\rho');
title('Zoom-in view near the contact wave');
subplot(2, 2, 3);
plot(xx, rhoEXACT, xx, rhoWENO, xx, rhoFWENO);
axis([0.15, 0.2, 0.1, 0.3]);
xlabel('x');
ylabel('\rho');
title('Zoom-in view near the shock wave');
subplot(2, 2, 4);
plot(0, 0, 0, 0, 0, 0);
legend('exact', 'WENO5', 'FWENO5', 'Location', 'southwest');
axis off;
print('-depsc', 'Sod2.eps');
