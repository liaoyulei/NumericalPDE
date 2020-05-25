%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate u_j^{n+1} for FWENO                             %
% u_j=u_j-\Delta t/\Delta x(\hat{f}_{j+1/2}-\hat{f}_{j-1/2}) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = FWENO (u, fu, alpha, dx, dt)
  epsilon = 1e-6;
  N = length(u)-6;
  p0 = zeros(1, N+6); %p_{3,0}(x_{j+1/2})
  p1 = zeros(1, N+6); %p_{3,1}(x_{j+1/2})
  p2 = zeros(1, N+6); %p_{3,2}{x_{j+1/2}}
  I0 = zeros(1, N+6); %I_{3,0}
  I1 = zeros(1, N+6); %I_{3,1}
  I2 = zeros(1, N+6); %I_{3,2}
  alpha0 = zeros(1, N+6); %\alpha_{3,0}
  alpha1 = zeros(1, N+6); %\alpha_{3,1}
  alpha2 = zeros(1, N+6); %\alpha_{3,2}
  omega0 = zeros(1, N+6); %\omega_{3,0}
  omega1 = zeros(1, N+6); %\omega_{3,1}
  omega2 = zeros(1, N+6); %\omega_{3,2}
  hatf = zeros(1, N+6);
  f = (fu + alpha .* u) ./ 2;
  p0(3: N+3)=(2 .* f(1: N+1) - 7 .* f(2: N+2) + 11 .* f(3: N+3)) ./ 6;
  p1(3: N+3) = (-f(2: N+2) + 5 .* f(3: N+3) + 2 .* f(4: N+4)) ./ 6;
  p2(3: N+3) = (2 .* f(3: N+3) + 5 .* f(4: N+4) - f(5: N+5)) ./ 6;
  I0(3: N+3) = (f(2: N+2) - f(3: N+3)).^2 + (f(3: N+3) - f(2: N+2)).^2;
  I1(3: N+3) = (f(3: N+3) - f(2: N+2)).^2 + (f(4: N+4) - f(3: N+3)).^2;
  I2(3: N+3) = (f(4: N+4) - f(3: N+3)).^2 + (f(5: N+5) - f(4: N+4)).^2;
  d(3: N+3) = (f(1: N+1) - 4 .* f(2: N+2) + 6 .* f(3: N+3) - 4 .* f(4: N+4) + f(5: N+5)).^2;
  alpha0(3: N+3) = 1/10 .* (1 + d(3: N+3) ./ (epsilon + I0(3: N+3)));
  alpha1(3: N+3) = 3/5 .* (1 + d(3: N+3) ./ (epsilon + I1(3: N+3)));
  alpha2(3: N+3) = 3/10 .* (1 + d(3: N+3) ./ (epsilon + I2(3: N+3)));
  omega0(3: N+3) = alpha0(3: N+3) ./ (alpha0(3: N+3) + alpha1(3: N+3) + alpha2(3: N+3));
  omega1(3: N+3) = alpha1(3: N+3) ./ (alpha0(3: N+3) + alpha1(3: N+3) + alpha2(3: N+3));
  omega2(3: N+3) = alpha2(3: N+3) ./ (alpha0(3: N+3) + alpha1(3: N+3) + alpha2(3: N+3)); 
  hatf(3: N+3) = omega0(3: N+3) .* p0(3: N+3) + omega1(3: N+3) .* p1(3: N+3) + omega2(3: N+3) .* p2(3: N+3); %\Hat{f}_{i+1/2},i=0,...N
  f = (fu - alpha .* u) ./ 2;
  p0(4: N+4)= (-f(2: N+2) + 5 .* f(3: N+3) + 2 .* f(4: N+4)) ./ 6;
  p1(4: N+4) = (2 .* f(3: N+3) + 5 .* f(4: N+4) - f(5: N+5)) ./ 6;
  p2(4: N+4) = (11 .* f(4: N+4) - 7 .* f(5: N+5) + 2 .* f(6: N+6)) ./ 6;
  I0(4: N+4) = (f(3: N+3) - f(2: N+2)).^2 + (f(4: N+4) - f(3: N+3)).^2;
  I1(4: N+4) = (f(4: N+4) - f(3: N+3)).^2 + (f(5: N+5) - f(4: N+4)).^2;
  I2(4: N+4) = (f(5: N+5) - f(4: N+4)).^2 + (f(6: N+6) - f(5: N+5)).^2;
  d(4: N+4) = (f(2: N+2) - 4 .* f(3: N+3) + 6 .* f(4: N+4) - 4 .* f(5: N+5) + f(6: N+6)).^2;
  alpha0(4: N+4) = 3/10 .* (1 + d(4: N+4) ./ (epsilon + I0(4: N+4)));
  alpha1(4: N+4) = 3/5 .* (1 + d(4: N+4) ./ (epsilon + I1(4: N+4)));
  alpha2(4: N+4) = 1/10 .* (1 + d(4: N+4) ./ (epsilon + I2(4: N+4)));
  omega0(4: N+4) = alpha0(4: N+4) ./ (alpha0(4: N+4) + alpha1(4: N+4) + alpha2(4: N+4));
  omega1(4: N+4) = alpha1(4: N+4) ./ (alpha0(4: N+4) + alpha1(4: N+4) + alpha2(4: N+4));
  omega2(4: N+4) = alpha2(4: N+4) ./ (alpha0(4: N+4) + alpha1(4: N+4) + alpha2(4: N+4)); 
  hatf(3: N+3) = hatf(3: N+3) + omega0(4: N+4) .* p0(4: N+4) + omega1(4: N+4) .* p1(4: N+4) + omega2(4: N+4) .* p2(4: N+4); %\Hat{f}_{i+1/2},i=0,...N
  u(4: N+3) = u(4: N+3) - dt./dx .* (hatf(4: N+3) - hatf(3: N+2));
  u(1: 3) = u(4: 6);
  u(N+4: N+6) = u(N+1: N+3);
end
