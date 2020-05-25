%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate the numerical solution %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = NuSolu (u0, xx, dx, dt, nt, format) %u0(1)=u0(nx)
  nx = length(xx); %length(u0)=length(xx)
  u = [u0, u0(2)]; %extend by periodic
  for it = 1: nt
    hatf = zeros(1, nx); %hatf=\hat{f}_{j+1/2},j=1:nx
    for j = 1: nx
      if strcmp(format, 'Lax-Friedrichs')
        hatf(j) = 1/2 * (1/2 * u(j)^2 + 1/2 * u(j+1)^2 - dx/dt * (u(j+1) - u(j)));
      elseif strcmp(format, 'Roe')
        if u(j) + u(j+1) > 0
          hatf(j) = 1/2 * u(j)^2;
        else
          hatf(j) = 1/2 * u(j+1)^2; 
        end
      elseif strcmp(format, 'Engquist-Osher')
        if u(j) > 0
          hatf(j) = hatf(j) + 1/2 * u(j)^2;
        end
        if u(j+1) < 0
          hatf(j) = hatf(j) + 1/2 * u(j+1)^2;
        end
      elseif strcmp(format, 'Godunov')
        if u(j) < u(j+1) && u(j) * u(j+1) > 0
          hatf(j) = 1/2 * min(u(j)^2, u(j+1)^2);
        elseif u(j) >= u(j+1)
          hatf(j) = 1/2 * max(u(j)^2, u(j+1)^2);
        end
      elseif strcmp(format, 'Lax-Wendroff')
        hatf(j) = 1/2 * (1/2 * u(j)^2 + 1/2 * u(j+1)^2 - dt/dx * 1/2 * (u(j) + u(j+1)) * (1/2 * u(j+1)^2 - 1/2 * u(j)^2));
      else
        break;
      end
    end
    u(2: nx) = u(2: nx) - dt / dx * (hatf(2: nx) - hatf(1: nx-1));
    u(1) = u(nx); %extend by periodic
    u(nx + 1) = u(2);
  end
  u = u(1: nx);
end
