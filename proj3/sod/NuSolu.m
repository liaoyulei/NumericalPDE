%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the numerical solution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho, v, p] = NuSolu (rho0, v0, p0, dx, dt, T, format)
  N = length(rho0); %length(u0)
  rho = [rho0(1: 3), rho0, rho0(N-2: N)];
  v = [v0(1: 3), v0, v0(N-2: N)];
  p = [p0(1: 3), p0, p0(N-2: N)];
  while T > 0
    dt = min(dt, T);
    [rho1, v1, p1] = comp_wise(rho, v, p, dx, dt, format);
    [rho1, v1, p1] = comp_wise(rho1, v1, p1, dx, dt, format);
    rho2 = 3/4 .* rho + 1/4 .* rho1;
    v2 = 3/4 .* v + 1/4 .* v1;
    p2 = 3/4 .* p + 1/4 .* p1;
    [rho2, v2, p2] = comp_wise(rho2, v2, p2, dx, dt, format);
    rho = 1/3 .* rho + 2/3 .* rho2;
    v = 1/3 .* v + 2/3 .* v2;
    p = 1/3 .* p + 2/3 .* p2;
    T = T - dt;
  end
  rho = rho(4: N+3);
  v = v(4: N+3);
  p = p(4: N+3);
end
