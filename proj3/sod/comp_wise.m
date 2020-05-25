%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Component-wise FD 1D system %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho, v, p] = comp_wise (rho, v, p, dx, dt, format)
  gamma = 1.4;
  str_eval = strcat(format, '(u, f, alpha)');
  rhov = rho .* v;
  e = p ./ (gamma - 1) + 1/2 .* rhov .* v;
  alpha = max(abs(v)) + max(gamma .* p ./ rho).^(1/2);
  rho = eval([format, '(rho, rhov, alpha, dx, dt)']);
  rhov = eval([format, '(rhov, rhov.*v+p, alpha, dx, dt)']);
  e = eval([format, '(e, v.*(e+p), alpha, dx, dt)']);
  v = rhov ./ rho;
  p = (gamma - 1) .* (e - 1/2 .* rhov .* v);
end
