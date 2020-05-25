%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the exact solution %
% f=u+sin(\pi(x-ut))=0         %
% by Newton's iteration method %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = ExSolu (x, t) %x array, t constant
  tgv = 1e-10; %error for Newton's iteration method
  nx = length(x);
  idx0 = (nx + 1) / 2;
  u = zeros(1, nx);
  u(1: idx0-1) = 1;
  u(idx0+1: nx) = -1;
  fu = u + sin(pi*(x-u*t))
  while max(abs(fu)) > tgv
    dfdu = 1 - pi * t * cos(pi*(x-u*t))
    u = u - fu ./ dfdu
    fu = u + sin(pi*(x-u*t))
  end
end
