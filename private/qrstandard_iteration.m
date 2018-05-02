function [H, G] = qrstandard_iteration(H, shift)
% QRSTANDARD_ITERATION    Perform a single step of QR iteration with shift
%
% This version of QR itertaion performs on an upper Hessenberg matrix 'H' with
% single-shifted algorithm.
%
% argin:
%   H     - An upper Hessenberg matrix to iterate
%   shift - A string of shift type (options: 'none', 'rayleigh' or 'wilkinson';
%           default:'wilkinson')
%
% argout:
%   H - Iteration result of H
%   G - The unitray matrices of the iteration H = P'HP, where
%       P = prod[G(:, :, k)]
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-03-31
% -------------------------------------------------

if nargin == 1
    shift = 'wilkinson';
end

n = length(H);

switch shift
    case 'none'
        mu = 0;
    case 'rayleigh'
        mu = H(n, n);
    case 'wilkinson'
        mu = qrstandard_wilkinson(H(n-1, n-1), H(n-1, n), H(n, n-1), H(n, n));
end

H(1, 1) = H(1, 1) - mu;
G = zeros(2, 2, n-1);

for k = 1:n-1
    [c, s, eta] = givens(H(k, k), H(k+1, k));
    G(:, :, k) = [c, s; -conj(s), conj(c)];
    H(k+1, k+1) = H(k+1, k+1) - mu;
    H(k:k+1, k) = [eta; 0];
    H(k:k+1, k+1:n) = G(:, :, k) * H(k:k+1, k+1:n);
end
for k = 1:n-1
    H(1:k+1, k:k+1) = H(1:k+1, k:k+1) * G(:, :, k)';
    H(k, k) = H(k, k) + mu;
end

H(n, n) = H(n, n) + mu;