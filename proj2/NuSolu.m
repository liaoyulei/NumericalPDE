%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the numerical solution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = NuSolu (u0, xx, dx, dt, nt, format, eq) %u0(1)=u0(nx)
  nx = length(xx); %length(u0)=length(xx)
  if eq == 1
    u = [u0(nx-2: nx-1), u0, u0(2: 3)]; %extend by periodic
  elseif eq == 2
    u = [u0(1), u0(1), u0, u0(nx), u0(nx)];
  end
  hatf = zeros(1, nx+4); %hatf(j)=\hatf{j+1/2}, j=2:nx+2
  u1 = zeros(1, nx+4);
  for it = 1: nt
    if format == 1 %Lax-Friedrichs
      hatf(2: nx+2) = 1/2 .* (1/2 .* u(2: nx+2).^2 + 1/2 * u(3: nx+3).^2 - (u(3: nx+3) - u(2: nx+2)));
      u(3: nx+2) = u(3: nx+2) - dt./dx .* (hatf(3: nx+2) - hatf(2: nx+1)); 
    elseif format == 2 %Godunov
      hatf(2: nx+2) = 1/2 .* (u(2: nx+2) < u(3: nx+3)) .* (u(2: nx+2) .* u(3: nx+3) > 0) .* min(u(2: nx+2).^2, u(3: nx+3).^2);
      hatf(2: nx+2) = hatf(2: nx+2) + 1/2 .* (u(2: nx+2) >= u(3: nx+3)) .* max(u(2: nx+2).^2, u(3: nx+3).^2);
      u(3: nx+2) = u(3: nx+2) - dt./dx .* (hatf(3: nx+2) - hatf(2: nx+1));
    elseif format == 3 %MUSCL
      u1 = MUSCL(u, dx, dt, eq);
      u = 1/2 .* u + 1/2 .* MUSCL(u1, dx, dt, eq);
    elseif format == 4 %Kurganov-Tadmor
      u1 = KT(u, dx, dt, eq);
      u1 = 3/4 .* u + 1/4 .* KT(u1, dx, dt, eq);
      u = 1/3 .* u + 2/3 .* KT(u1, dx, dt, eq);
    end
    if eq == 1
      u(1: 2) = u(nx: nx+1); %extend by periodic
      u(nx+3: nx+4) = u(4: 5);
    end
  end
  u = u(3: nx+2);
end
