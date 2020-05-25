%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the numerical solution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = NuSolu (u0, dx, dt, T, format, eq, init)
  N = length(u0); %length(u0)
  if init == 4
    u = [u0(1: 3), u0, u0(N-2: N)];
  else
    u = [u0(N-2: N), u0, u0(1: 3)]; %extend by periodic
  end
  if strcmp(format, 'WENO5')
    while T > dt
      u1 = WENO(u, dx, dt, eq, init);
      u1 = 3/4 .* u + 1/4 .* WENO(u1, dx, dt, eq, init);
      u = 1/3 .* u + 2/3 .* WENO(u1, dx, dt, eq, init);
      T = T - dt;
    end
    u1 = WENO(u, dx, T, eq, init);
    u1 = 3/4 .* u + 1/4 .* WENO(u1, dx, T, eq, init);
    u = 1/3 .* u + 2/3 .* WENO(u1, dx, T, eq, init);
  elseif strcmp(format, 'FWENO5')
    while T > dt
      u1 = FWENO(u, dx, dt, eq, init);
      u1 = 3/4 .* u + 1/4 .* FWENO(u1, dx, dt, eq, init);
      u = 1/3 .* u + 2/3 .* FWENO(u1, dx, dt, eq, init);
      T = T - dt;
    end
    u1 = FWENO(u, dx, T, eq, init);
    u1 = 3/4 .* u + 1/4 .* FWENO(u1, dx, T, eq, init);
    u = 1/3 .* u + 2/3 .* FWENO(u1, dx, T, eq, init);
  else
    if strcmp(format, 'ESWENO5')
      phi = 0;
    elseif strcmp(format, 'ESWENO6')
      phi = 1/20;
    end
    while T > dt
      u1 = ESWENO(u, dx, dt, eq, init, phi);
      u1 = 3/4 .* u + 1/4 .* ESWENO(u1, dx, dt, eq, init, phi);
      u = 1/3 .* u + 2/3 .* ESWENO(u1, dx, dt, eq, init, phi);
      T = T - dt;
    end
    u1 = ESWENO(u, dx, T, eq, init, phi);
    u1 = 3/4 .* u + 1/4 .* ESWENO(u1, dx, T, eq, init, phi);
    u = 1/3 .* u + 2/3 .* ESWENO(u1, dx, T, eq, init, phi);
  end
  u = u(4: N+3);
end
