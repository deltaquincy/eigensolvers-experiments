function mu = qrstandard_wilkinson(alpha, beta, gamma, delta)
% QRSTANDARD_WILKINSON    Compute Wilkinson's shift of standard QR algorithm
%
% For a 2x2 matrix [alpha, beta; gamma, delta], Wilkinson's shift means the
% eigenvalue of this matrix which is closer to 'delta'.
%
% argin:
%   alpha, beta, gamma, delta - Elements of the 2x2 matrix to compute the shift
%
% argout:
%   mu - Wilkinson's shift of this matrix
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-03-31
% -------------------------------------------------

mu = delta;
s = sum(abs([alpha, beta, gamma, delta]));

if s == 0
    return
end

q = (beta / s) * (gamma / s);

if q ~= 0
    p = 0.5 * (alpha / s - delta / s);
    r = sqrt(p^2 + q);
    if real(p)*real(r) + imag(p)*imag(r) < 0
        r = -r;
    end
    mu = mu - s*(q/(p+r));
end

