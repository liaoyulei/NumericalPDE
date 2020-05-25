%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eq='Advection': u_t+u_x=0 -1\le x\le 1    %
% eq='Burgers': u_t+(u^2/2)_x=0 -1\le x\le 1%
% u_1(x,0)=\sin(\pi x)                      %
% u_2(x,0)=\sin^4(\pi x)                    %
% u_3(x,0)=0.25+0.5\sin(\pi x)              %
% u_4(x,0)=-1 x<0, 1x>0                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = ExSolu (x, t, eq, init) %x array, t constant, eq equation
  eps = 1e-8;
  if strcmp(eq, 'Advection')
    if init == 1
      u = sin(pi.*(x-t));
    elseif init == 2
      u = sin(pi.*(x-t)).^4;
    elseif init == 3
      u = 1/4 + 1/2 .* sin(pi.*(x-t));
    end
  elseif strcmp(eq, 'Burgers')
    if init == 3
      if t < 10
        u = 1/4 + 1/2 .* sin(pi.*x);
      else
        u = 0.08 .* x + (0.25 - 0.08 .* sign(x));
      end
      fu = u - 1/4 - 1/2 .* sin(pi.*(x-u.*t));
      while max(abs(fu)) > eps
        dfdu = 1 + 1/2 .* pi .* t .* cos(pi.*(x-u.*t));
        u = u - fu ./ dfdu;
        fu = u - 1/4 - 1/2 .* sin(pi.*(x-u.*t));
      end
    elseif init == 4
      if t == 0
        u = zeros(1, length(x));
        u(x<0) = -1;
        u(x>0) = 1;
      else
        u = x ./ t;
        u(x<-t) = -1;
        u(x>t) = 1;
      end
    end
  end
end
