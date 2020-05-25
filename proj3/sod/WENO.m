%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate u_j^{n+1} for WENO                               %
% u_j=u_j-\Delta t/\Delta x(\hat{f}_{j+1/2}-\hat{f}_{j-1/2}) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = WENO (u, fu, alpha, dx, dt)
  epsilon = 1e-6;
  N = length(u)-6;
  v0 = zeros(1, N+6);
  v1 = zeros(1, N+6);
  v2 = zeros(1, N+6);
  beta0 = zeros(1, N+6);
  beta1 = zeros(1, N+6);
  beta2 = zeros(1, N+6);
  alpha0 = zeros(1, N+6);
  alpha1 = zeros(1, N+6);
  alpha2 = zeros(1, N+6);
  omega0 = zeros(1, N+6);
  omega1 = zeros(1, N+6);
  omega2 = zeros(1, N+6);
  hatf = zeros(1, N+6);
  v = (fu + alpha .* u) ./ 2;
  v0(3: N+4) = 1/3 .* v(3: N+4) + 5/6 .* v(4: N+5) - 1/6 .* v(5: N+6); %v_{i+1/2}^{(0)}
  v1(3: N+4) = -1/6 .* v(2: N+3) + 5/6 .* v(3: N+4) + 1/3 .* v(4: N+5); %v_{i+1/2}^{(1)}
  v2(3: N+4) = 1/3 .* v(1: N+2) - 7/6 .* v(2: N+3) + 11/6 .* v(3: N+4); %v_{i+1/2}^{(2)}
  beta0(3: N+4) = 13/12 .* (v(3: N+4) - 2 .* v(4: N+5) + v(5: N+6)).^2 + 1/4 .* (3 .* v(3: N+4) - 4 .* v(4: N+5) + v(5: N+6)).^2; %\beta_0
  beta1(3: N+4) = 13/12 .* (v(2: N+3) - 2 .* v(3: N+4) + v(4: N+5)).^2 + 1/4 .* (v(2: N+3) - v(4: N+5)).^2; %\beta_1
  beta2(3: N+4) = 13/12 .* (v(1: N+2) - 2 .* v(2: N+3) + v(3: N+4)).^2 + 1/4 .* (v(1: N+2) - 4 .* v(2: N+3) + 3 .* v(3: N+4)).^2; %\beta_2
  alpha0(3: N+4) = 3/10 ./ (epsilon + beta0(3: N+4)).^2; %\alpha_r=d_r/(\epsilon+\beta_r)^2
  alpha1(3: N+4) = 3/5 ./ (epsilon + beta1(3: N+4)).^2;
  alpha2(3: N+4) = 1/10 ./ (epsilon + beta2(3: N+4)).^2;  
  omega0(3: N+4) = alpha0(3: N+4) ./ (alpha0(3: N+4) + alpha1(3: N+4) + alpha2(3: N+4)); %\omega_r=\alpha_r/(\sum_{s=0}^2\alpha_s)
  omega1(3: N+4) = alpha1(3: N+4) ./ (alpha0(3: N+4) + alpha1(3: N+4) + alpha2(3: N+4));
  omega2(3: N+4) = alpha2(3: N+4) ./ (alpha0(3: N+4) + alpha1(3: N+4) + alpha2(3: N+4));
  hatf(3: N+3) = omega0(3: N+3) .* v0(3: N+3) + omega1(3: N+3) .* v1(3: N+3) + omega2(3: N+3) .* v2(3: N+3); %v_{i+1/2}^{-}
  v = (fu - alpha .* u) ./ 2;
  v0(3: N+4) = 11/6 .* v(3: N+4) - 7/6 .* v(4: N+5) + 1/3 .* v(5: N+6); %u_{i-1/2}^{(0)}=11/6u_i-7/6u_{i+1}+1/3u_{i+1}
  v1(3: N+4) = 1/3 .* v(2: N+3) + 5/6 .* v(3: N+4) - 1/6 .* v(4: N+5); %u_{i-1/2}^{(1)}=1/3u_{i-1}+5/6u_i-1/6u_{i+1}
  v2(3: N+4) = -1/6 .* v(1: N+2) + 5/6 .* v(2: N+3) + 1/3 .* v(3: N+4); %u_{i-1/2}^{(2)}=-1/6u_{i-2}+5/6u_{i-1}+1/3u_i
  beta0(3: N+4) = 13/12 .* (v(3: N+4) - 2 .* v(4: N+5) + v(5: N+6)).^2 + 1/4 .* (3 .* v(3: N+4) - 4 .* v(4: N+5) + v(5: N+6)).^2; %\beta_0
  beta1(3: N+4) = 13/12 .* (v(2: N+3) - 2 .* v(3: N+4) + v(4: N+5)).^2 + 1/4 .* (v(2: N+3) - v(4: N+5)).^2; %\beta_1
  beta2(3: N+4) = 13/12 .* (v(1: N+2) - 2 .* v(2: N+3) + v(3: N+4)).^2 + 1/4 .* (v(1: N+2) - 4 .* v(2: N+3) + 3 .* v(3: N+4)).^2; %\beta_2
  alpha0(3: N+4) = 1/10 ./ (epsilon + beta0(3: N+4)).^2; %\Tilde{\alpha}_r=\Tilde{d}_r/(\epsilon+\beta_r)^2
  alpha1(3: N+4) = 3/5 ./ (epsilon + beta1(3: N+4)).^2;
  alpha2(3: N+4) = 3/10 ./ (epsilon + beta2(3: N+4)).^2;
  omega0(3: N+4) = alpha0(3: N+4) ./ (alpha0(3: N+4) + alpha1(3: N+4) + alpha2(3: N+4)); %\omega_r=\alpha_r/(\sum_{s=0}^2\alpha_s)
  omega1(3: N+4) = alpha1(3: N+4) ./ (alpha0(3: N+4) + alpha1(3: N+4) + alpha2(3: N+4));
  omega2(3: N+4) = alpha2(3: N+4) ./ (alpha0(3: N+4) + alpha1(3: N+4) + alpha2(3: N+4));
  hatf(3: N+3) = hatf(3: N+3) + omega0(4: N+4) .* v0(4: N+4) + omega1(4: N+4) .* v1(4: N+4) + omega2(4: N+4) .* v2(4: N+4); %u_{i+1/2}^{+}
  u(4: N+3) = u(4: N+3) - dt./dx .* (hatf(4: N+3) - hatf(3: N+2));
  u(1: 3) = u(4: 6);
  u(N+4: N+6) = u(N+1: N+3);
end
