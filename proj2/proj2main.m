%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve Burgers equation         %
% u_t+(u^2/2)_x=0                %
% u_1(x,0)=sin(\pi x), x\in[0,2] %
% u_2(x,0)=1 x\le0,\\-0.5 x>0    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%choose one of the five formats and the equation
%format = 1; %Lax-Friedrichs;
%format = 2; %Godunov;
format = 3; %MUSCL;
%format = 4; %Kurganov-Tadmor;
string = {'Lax-Friedrichs', 'Godunov', 'MUSCL', 'Kurganov-Tadmor'};
%eq = 1;
eq = 2;
t0 = 0; %init time
if eq == 1
  x1 = 0; %init space
  x2 = 2; %final space
elseif eq == 2
  x1 = -1; %init space
  x2 = 1; %final space
end
r = 1/4; %dt/dx

% task1 calculate error and order
ne = 6; %number of test
e = zeros(1, ne);
order = zeros(1, ne);
tf = 0.15; %final time
nt = 6; % number of time interval in the first
for in = 1: ne
  nt = 2 .* nt;
  dt = (tf - t0) ./ nt;
  dx = dt ./ r;
  xx = x1: dx: x2; %make sure xx[nx]=x2
  u0= ExSolu(xx, 0, eq);
  uex = ExSolu(xx, tf, eq);
  u = NuSolu(u0, xx, dx, dt, nt, format, eq);
  e(in) = sum(abs(u-uex).*dx);
  if in > 1
    order(in) = log2(e(in-1)./e(in));
  end
end

%task2 draw exact and numerical solution
tf = 0.5; %final time
nt = 2^8;
dt = (tf - t0) ./ nt;
dx = dt ./ r;
xx = x1: dx: x2;
u0 = ExSolu(xx, 0, eq); %make sure xx[nx]=x2
uex = ExSolu(xx, tf, eq);
u = NuSolu(u0, xx, dx, dt, nt, format, eq);
%subplot(2, 2, format);
plot(xx, uex, xx, u);
axis([x1, x2, -1.1, 1.1]);
xlabel('x');
ylabel('u');
gca = legend('exact', 'numerical', 'Location','southwest');
%po = get(gca, 'Position');
%set(gca, 'Position', [po(1)-0.03, po(2), po(3), po(4)]);
title(strcat(string(format), ' scheme: t=0.5, r=1/4'));
%print('-depsc', 'name.eps');
