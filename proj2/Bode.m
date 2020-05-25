
function u = Bode (x, t, dx, eq) %x array, t constant, eq equation
  nx = length(x);
  u = zeros(1, nx);
  u = ExSolu(x, t, eq);
%  u(2: nx-1) = 24 .* u(2: nx-1);
%  u(2: nx-1) = u(2: nx-1) + 14 .* ExSolu(x(2: nx-1)-1/2*dx, t, eq);
%  u(2: nx-1) = u(2: nx-1) + 64 .* ExSolu(x(2: nx-1)-1/4*dx, t, eq);
%  u(2: nx-1) = u(2: nx-1) + 64 .* ExSolu(x(2: nx-1)+1/4*dx, t, eq);
%  u(2: nx-1) = u(2: nx-1) + 14 .* ExSolu(x(2: nx-1)+1/2*dx, t, eq);
%  u(2: nx-1) = u(2: nx-1) ./ 180;
  for j = 1: 100
    u(2: nx-1) = u(2: nx-1) + ExSolu(x(2: nx-1)-j/200*dx, t, eq) + ExSolu(x(2: nx-1)+j/200*dx, t,eq);
  end
  u(2: nx-1) = u(2: nx-1) / 201;
end
