

function [circAvg, vecNorm, vecAngle, unitVector] = circularAverage(F, angles, cosang, sinang)
% CIRCULARAVERAGE calculates the circular average
%       [circAvg, vecNorm, vecAngle, unitVector] = circularAverage( F, angles )
%
%       Z = Fexp(it) = F*cos(t) + i*F*sin(t);
%   circAvg = <Z>_t

if nargin < 4
    cosang = cos(angles);
    sinang = sin(angles);
end

if max(F) >= 1
    %      ReZ = F .* cosang;
    %      ImZ = F .* sinang;
    
    sF = sum(F);
    circAvg = [sum(F .* cosang), sum(F .* sinang)]/sF;
else
    circAvg = [NaN NaN];
end

if nargout > 1
    vecNorm = norm(circAvg);
    vecAngle = atan2( circAvg(2), circAvg(1) );
    unitVector = circAvg ./ vecNorm;
end

end
