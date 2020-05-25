%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the exact solution %
% For eq1                      %
% f=u-sin(\pi(x-ut))=0         %
% by Newton's iteration method %
% For eq2                      %
% u=1,x\le t/4\\-0.5,x>t/4     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = ExSolu (x, t, eq) %x array, t constant, eq equation
  nx = length(x);
  u = zeros(1, nx);
  if eq == 1
    tgv = 1e-10; %error for Newton's iteration method 
    if t < 0.5
      u = sin(pi.*x);
    else
      u = sign(1-x) .* sin(pi.*x./2);
    end
    fu = u - sin(pi.*(x-u.*t));
    while max(abs(fu)) > tgv
      dfdu = 1 + pi .* t .* cos(pi.*(x-u.*t));
      u = u - fu ./ dfdu;
      fu = u - sin(pi.*(x-u.*t));
    end
  elseif eq == 3
    u(x<=t/4) = 1;
    u(x>t/4) = -0.5;
  elseif eq == 2
    if t == 0
      u(x<0) = -1;
      u(x>=0) = 1;
    else
      u = x ./ t;
      u(x<-t) = -1;
      u(x>t) = 1;
  end
end
