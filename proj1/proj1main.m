%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve Burgers equation           %
% u_t+(u^2/2)_x=0, x\in[-1,1], t>0 %
% u(x, 0)=sin(\pi x+\pi) periodic  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose one of the five formats
%format = 'Lax-Friedrichs';
%format = 'Roe';
%format = 'Engquist-Osher';
%format = 'Godunov';
format = 'Lax-Wendroff';
t0 = 0; %init time
x1 = -1; %init space
x2 = 1; %final space
r = 1; %dt/dx

% task1 calculate error and order
ne = 6; %number of test
e = zeros(1, ne);
order = zeros(1, ne);
tf = 0.15; %final time
nt = 3; % number of time interval in the first
for in = 1: ne
  nt = nt * 2;
  dt = (tf - t0) / nt;
  dx = dt / r;
  xx = x1: dx: x2; %make sure xx[nx]=x2
  u0 = -sin(pi*xx);
  uex = ExSolu(xx, tf);
  u = NuSolu(u0, xx, dx, dt, nt, format);
  e(in) = sum(abs(u-uex)*dx);
  if in > 1
    order(in) = log2(e(in-1)/e(in));
  endif
end

%task2 draw exact and numerical solution
tf = 0.5; %final time
nt = 2^7;
dt = (tf - t0) / nt;
dx = dt / r;
xx = x1: dx: x2;
u0 = -sin(pi*xx); %make sure xx[nx]=x2
uex = ExSolu(xx, tf);
u = NuSolu(u0, xx, dx, dt, nt, format);
plot(xx, uex, xx, u);
axis([x1, x2, -1.1, 1.1]);
xlabel('x');
ylabel('u');
legend('exact', 'numerical');
title(strcat(format, ' scheme: t=0.5, r=1'));
print('-depsc', strcat(format, '.eps'));
