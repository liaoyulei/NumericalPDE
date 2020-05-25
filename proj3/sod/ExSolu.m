%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exact solution for u_t+f(u)_x=0                   %
% u=(\rho,\rho v,E)^T,f=(\rho v,\rho v^2+p,(E+p))^T %
% (\rho_0,v_0,p_0)=(1,0,1) x<0=\\(0.125,0,0.1) x>0  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rho, v, p] = ExSolu (xx, t);
  N = length(xx);
  rho = zeros(1, N);
  v = zeros(1, N);
  p = zeros(1, N);
  if t == 0
    rho(xx<=0) = 1;
    rho(xx>0) = 0.125;
    p(xx<=0) = 1;
    p(xx>0) = 0.1;
  else
    v = 5/6 .* (xx ./ t + 1.4^(1/2));
    c = v - xx ./ t;
    p = (c ./ 1.4^(1/2)).^7;
    rho = (c ./ 1.4^(1/2)).^5;
    rho(xx<=-1.1832159566199*t) = 1;
    v(xx<=-1.1832159566199*t) = 0;
    p(xx<=-1.1832159566199*t) = 1;
    rho(xx>-0.0702728125612*t) = 0.4263194281785;
    v(xx>-0.0702728125612*t) = 0.9274526200490;
    p(xx>-0.0702728125612*t) = 0.3031301780506;
    rho(xx>0.9274526200490*t) = 0.2655737117053;
    v(xx>0.9274526200490*t) = 0.9274526200490;
    p(xx>0.9274526200490*t) = 0.3031301780506;
    rho(xx>1.7521557320302*t) = 0.125;
    v(xx>1.7521557320302*t) = 0;
    p(xx>1.7521557320302*t) = 0.1;
  end
end
