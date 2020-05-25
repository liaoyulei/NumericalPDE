%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate u_j^{n+1} for ESWENO                             %
% u_j=u_j-\Delta t/\Delta x(\hat{f}_{j+1/2}-\hat{f}_{j-1/2}) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = ESWENO (u, dx, dt, eq, init, phi)
  N = length(u)-6;
  fLL = zeros(1, N+6); %f^{LL}(u_{j+1/2})
  fL = zeros(1, N+6); %f^L(u_{j+1/2})
  fR = zeros(1, N+6); %f^R(u_{j+1/2})
  fRR = zeros(1, N+6); %f^{RR}(u_{j+1/2})
  betaLL = zeros(1, N+6); %\beta_{LL}
  betaL = zeros(1, N+6); %\beta_L
  betaR = zeros(1, N+6); %\beta_R
  betaRR = zeros(1, N+6); %\beta_{RR}
  tau = zeros(1, N+6); %\tau_5
  alphaLL = zeros(1, N+6); %\alpha_{LL}
  alphaL = zeros(1, N+6); %\alpha_L
  alphaR = zeros(1, N+6); %\alpha_R
  alphaRR = zeros(1, N+6); %\alpha_{RR}
  wLL = zeros(1, N+6); %w_{j+1/2}^{(LL)}
  wL = zeros(1, N+6); %w_{j+1/2}^{(L)}
  wR = zeros(1, N+6); %w_{j+1/2}^{(R)}
  wRR = zeros(1, N+6); %w_{j+1/2}^{(RR)}
  hatf = zeros(1, N+6);
  if strcmp(eq, 'Advection')
    f = u;
  elseif strcmp(eq, 'Burgers')
    f = u.^2 ./ 2;
  end
  fLL(3: N+3)=(2 .* f(1: N+1) - 7 .* f(2: N+2) + 11 .* f(3: N+3)) ./ 6;
  fL(3: N+3) = (-f(2: N+2) + 5 .* f(3: N+3) + 2 .* f(4: N+4)) ./ 6;
  fR(3: N+3) = (2 .* f(3: N+3) + 5 .* f(4: N+4) - f(5: N+5)) ./ 6;
  fRR(3: N+3) = (11 .* f(4: N+4) - 7 .* f(5: N+5) + 2 .* f(6: N+6)) ./ 6;
  betaLL(3: N+3) = 13/12 .* (f(1: N+1) - 2 .* f(2: N+2) + f(3: N+3)).^2 + 1/4 .* (f(1: N+1) - 4 .* f(2: N+2) + 3 .* f(3: N+3)).^2;
  betaL(3: N+3) = 13/12 .* (f(2: N+2) - 2 .* f(3: N+3) + f(4: N+4)).^2 + 1/4 .* (f(2: N+2) - f(4: N+4)).^2;
  betaR(3: N+3) = 13/12 .* (f(3: N+3) - 2 .* f(4: N+4) + f(5: N+5)).^2 + 1/4 .* (3 .* f(3: N+3) - 4 .* f(4: N+4) + f(5: N+5)).^2;
  betaRR(3: N+3) = 13/12 .* (f(4: N+4) - 2 .* f(5: N+5) + f(6: N+6)).^2 + 1/4 .* (-5 .* f(4: N+4) + 8 .* f(5: N+5) - 3 .* f(6: N+6)).^2;
  if phi ~= 0
    tau(3: N+3) = (-f(1: N+1) + 5 .* f(2: N+2) - 10 .* f(3: N+3) + 10 .* f(4: N+4) - 5 .* f(5: N+5) + f(6: N+6)).^2;
    epsilon = dx.^6;
  else
    tau(3: N+3) = (f(1: N+1) - 4 .* f(2: N+2) + 6 .* f(3: N+3) - 4 .* f(4: N+4) + f(5: N+5)).^2;
    epsilon = dx.^5;
  end
  alphaLL(3: N+3) = (1/10 - phi) .* (1 + tau(3: N+3) ./ (epsilon + betaLL(3: N+3)));
  alphaL(3: N+3) = (6/10 - 3 .* phi) .* (1 + tau(3: N+3) ./ (epsilon + betaL(3: N+3)));
  alphaR(3: N+3) = (3/10 + 3 .* phi) .* (1 + tau(3: N+3) ./ (epsilon + betaR(3: N+3)));
  alphaRR(3: N+3) = phi .* (1 + tau(3: N+3) ./ (epsilon + betaRR(3: N+3)));
  wLL(3: N+3) = alphaLL(3: N+3) ./ (alphaLL(3: N+3) + alphaL(3: N+3) + alphaR(3: N+3) + alphaRR(3: N+3));
  wL(3: N+3) = alphaL(3: N+3) ./ (alphaLL(3: N+3) + alphaL(3: N+3) + alphaR(3: N+3) + alphaRR(3: N+3));
  wR(3: N+3) = alphaR(3: N+3) ./ (alphaLL(3: N+3) + alphaL(3: N+3) + alphaR(3: N+3) + alphaRR(3: N+3));
  wRR(3: N+3) = alphaRR(3: N+3) ./ (alphaLL(3: N+3) + alphaL(3: N+3) + alphaR(3: N+3) + alphaRR(3: N+3));  
  hatf(3: N+3) = wLL(3: N+3) .* fLL(3: N+3) + wL(3: N+3) .* fL(3: N+3) + wR(3: N+3) .* fR(3: N+3) + wRR(3: N+3) .* fRR(3: N+3); %\Hat{f}_{i+1/2},i=0,...N
  u(4: N+3) = u(4: N+3) - dt./dx .* (hatf(4: N+3) - hatf(3: N+2));
  if init == 4
    u(1: 3) = u(4: 6);
    u(N+4: N+6) = u(N+1: N+3);
  else
    u(1: 3) = u(N+1: N+3);
    u(N+4: N+6) = u(4: 6);
  end
end
