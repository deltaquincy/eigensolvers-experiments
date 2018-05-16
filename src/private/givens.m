function [c, s, eta] = givens(alpha, beta)
% GIVENS    Compute Givens rotation.
%
% The Givens rotation is an unitray transformation, from [a; b] to [eta; 0].
% This function computes the sine, cosine value and the scalar eta.
%
% argin:
%   alpha, beta - Two scalars in the vector to rotate (real or complex).
%
% argout:
%   usage 1: [c, s, eta] = givens(alpha, beta)
%       c, s - The cosine and sine values (scalars) in the Givens rotation.
%       eta  - The eta value (scalar) in the rotated result.
%   usage 2: G = givens(alpha, beta)
%       G - The 2*2 Givens rotation matrix that reduce [a; b] to [eta; 0].
%   
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-03-31
% -------------------------------------------------

if beta == 0
    c = 1;
    s = 0;
    eta = alpha;
    return
end

if alpha == 0
    c = 0;
    s = 1;
    eta = beta;
    return
end

abs_alpha = abs(alpha);  % Prevent to compute for three times.
mu = alpha / abs_alpha;
tau = abs_alpha + abs(beta);
delta = tau * sqrt(abs(alpha/tau)^2 + abs(beta/tau)^2);

c = abs_alpha / delta;
s = mu * conj(beta) / delta;
eta = delta * mu;

% Compute G if only one output is needed.
if nargout == 1
    c = [c, s; -conj(s), conj(c)];
end



