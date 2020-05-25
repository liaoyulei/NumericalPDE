%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minmod(ul,ur)=                  %
% a, \abs{ul}<\abs{ur} && ul*ur>0 %
% b, \abs{ul}<\abs{ur} && ul*ur>0 %
% 0, ul*ur\le0                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tildeu = minmod (ul, ur)
  tildeu = 1/2 .* (sign(ul) + sign(ur)) .* min(abs(ul), abs(ur));
end
