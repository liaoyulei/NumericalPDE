%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate u_j^{n+1} for extend MUSCL                       %
% u_j=u_j-\Delta t/\Delta x(\hat{f}_{j+1/2}-\hat{f}_{j-1/2}) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = KT (u, dx, dt, eq)
  nx = length(u)-4;
  uminus = zeros(1, nx+4); %uminus(j)=u^{-}_{j+1/2}, j=2:nx+3
  uplus = zeros(1, nx+4);  %uplus(j)=u^{+}_{j+1/2}, j=1:nx+2 
  umod1 = zeros(1, nx+4); %umod1(j)=\tilde{u}^{mod}_j, j=2:nx+3
  umod2 = zeros(1, nx+4); %umod2(j)=\tilde{\tilde{u}}^{mod}_j, j=2:nx+3
  hatf = zeros(1, nx+4); %hatf(j)=\hatf{j+1/2}, j=2:nx+2
  uminus(2: nx+3) = -1/6 .* u(1: nx+2) + 5/6 .* u(2: nx+3) + 1/3 .* u(3: nx+4);
  uplus(1: nx+2) = 1/3 .* u(1: nx+2) + 5/6 .* u(2: nx+3) -1/6 .* u(3: nx+4); 
  umod1(2: nx+3) = minmod(uminus(2: nx+3)-u(2: nx+3), minmod(u(3: nx+4)-u(2: nx+3), u(2: nx+3)-u(1: nx+2)));
  umod2(2: nx+3) = minmod(u(2: nx+3)-uplus(1: nx+2), minmod(u(3: nx+4)-u(2: nx+3), u(2: nx+3)-u(1: nx+2)));
  uminus(2: nx+3) = u(2: nx+3) + umod1(2: nx+3);
  uplus(1: nx+2) = u(2: nx+3) - umod2(2: nx+3);
%  hatf(2: nx+2) = 1/2 .* (1/2 .* uminus(2: nx+2).^2 + 1/2 * uplus(2: nx+2).^2 - (uplus(2: nx+2) - uminus(2: nx+2)));
  hatf(2: nx+2) = 1/2 .* (uminus(2: nx+2) < uplus(2: nx+2)) .* (uminus(2: nx+2) .* uplus(2: nx+2) > 0) .* min(uminus(2: nx+2).^2, uplus(2: nx+2).^2);
  hatf(2: nx+2) = hatf(2: nx+2) + 1/2 .* (uminus(2: nx+2) >= uplus(2: nx+2)) .* max(uminus(2: nx+2).^2, uplus(2: nx+2).^2);
  u(3: nx+2) = u(3: nx+2) - dt./dx .* (hatf(3: nx+2) - hatf(2: nx+1));
  if eq == 1
    u(1: 2) = u(nx: nx+1);
    u(nx+3: nx+4) = u(4: 5);
  end
end
