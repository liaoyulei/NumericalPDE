%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve conservation law equations %
% u_t+u_x=0 -1\le x\le 1           %
% u_1(x,0)=\sin(\pi x)             %
% u_2(x,0)=\sin^4(\pi x)           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%choose one format
%format = 'WENO5';
format = 'FWENO5';
%choose the equations
eq = 'Advection'; %choose init=1,2,3
%eq = 'Burgers'; % choose init=3,4
%choose the init value
%init = 1; %\sin(\pi x)
init = 2; %\sin^4(\pi x)
%init = 3; %0.25+0.5\sin(\pi x)
%init = 4; %-1 x<0,1 x>0
x1 = -1; %init space
x2 = 1; %final space
CFL = 1; %dt^3/dx^5

% task1 calculate error and order
ne = 6; %number of test
e1 = zeros(1, ne);
order1 = zeros(1, ne);
einfty = zeros(1, ne);
orderinfty = zeros(1, ne);
if strcmp(eq, 'Advection')
  T = 1; %final time
  N = 5; %number of space interval
elseif strcmp(eq, 'Burgers')
  T = 0.3;
  N = 5;
end
for in = 1: ne
  N = 2 .* N;
  dx = (x2 - x1) ./ N;
  dt = (CFL .* dx.^5).^(1/3);
  xx = x1+dx./2: dx: x2;
  u0= ExSolu(xx, 0, eq, init);
  uex = ExSolu(xx, T, eq, init);
  u = NuSolu(u0, dx, dt, T, format, eq, init);
  e1(in) = sum(abs(u-uex).*dx);
  einfty(in) = max(abs(u-uex));
  if in > 1
    order1(in) = log2(e1(in-1)./e1(in));
    orderinfty(in) = log2(einfty(in-1)./einfty(in));
  end
end

% draw solutions for Burgers equations
eq = 'Burgers';
init = 3;
x1 = -1; %init space
x2 = 1; %final space
CFL = 1; %dt^3/dx^5
T = 12;
N = 200; %number of space interval
dx = (x2 - x1) ./ N;
dt = (CFL .* dx.^5).^(1/3);
xx = x1+dx./2: dx: x2;
u0 = ExSolu(xx, 0, eq, init);
uEXACT = ExSolu(xx, T, eq, init);
uWENO = NuSolu(u0, dx, dt, T, 'WENO5', eq, init);
uFWENO =NuSolu(u0, dx, dt, T, 'FWENO5', eq, init);
subplot(1, 2, 1);
plot(xx, uEXACT, xx, uWENO, xx, uFWENO);
axis([-1, 1, 0.16, 0.34]);
xlabel('x');
ylabel('u');
title('T=12, N=200');
gca = legend('exact', 'WENO5', 'FWENO5');
po = get(gca, 'Position');
set(gca, 'Position', [po(1)+0.05, po(2), po(3), po(4)]);
subplot(1, 2, 2);
plot(xx, uEXACT, xx, uWENO, xx, uFWENO);
axis([-0.1, 0.4, 0.16, 0.22]);
xlabel('x');
ylabel('u');
title('Zoom-in view');
legend('exact', 'WENO5', 'FWENO5');
print('-depsc', 'Burgers1.eps');

% draw solutions for Burgers equations
eq = 'Burgers';
init = 4;
x1 = -1; %init space
x2 = 1; %final space
CFL = 1; %dt^3/dx^5
T = 0.8;
N = 200; %number of space interval
dx = (x2 - x1) ./ N;
dt = (CFL .* dx.^5).^(1/3);
xx = x1+dx./2: dx: x2;
u0 = ExSolu(xx, 0, eq, init);
uEXACT = ExSolu(xx, T, eq, init);
uWENO = NuSolu(u0, dx, dt, T, 'WENO5', eq, init);
uFWENO =NuSolu(u0, dx, dt, T, 'FWENO5', eq, init);
subplot(1, 2, 1);
plot(xx, uEXACT, xx, uWENO, xx, uFWENO);
axis([-1, 1, -1.1, 1.1]);
xlabel('x');
ylabel('u');
title('T=0.8, N=200');
legend('exact', 'WENO5', 'FWENO5', 'Location','northwest');
subplot(1, 2, 2);
plot(xx, uEXACT, xx, uWENO, xx, uFWENO);
axis([-1, -0.5, -1.1, -0.6]);
xlabel('x');
ylabel('u');
title('Zoom-in view');
legend('exact', 'WENO5', 'FWENO5', 'Location','northwest');
print('-depsc', 'Burgers2.eps');
