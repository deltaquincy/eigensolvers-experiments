function [H, U] = qrstandard_reduction(A)
% QRSTANDARD_REDUCTION    Compute upper Hessenberg reduction
%
% This function computes the upper Hessenberg reduction of a matrix A, so that
% U'AU = H. Note that A should be a square matrix.
%
% argin:
%   A - A matrix to perform upper Hessenberg reduction (real or complex)
%
% argout:
%   H - The upper Hessenberg matrix which is unitray similar with A
%   U - the unitray matrix transforming A to H, so that U'AU = H
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-03-31
% -------------------------------------------------

n = size(A);
if n(1) ~= n(2)
    error('Fatal Error: The input is not a square matrix!');
end
n = n(1);

H = A;
if nargout == 2
    U = zeros(n, n);
end

for k = 1:n-2
    w = house(H(k+1:n, k));
    v = w' * H(k+1:n, k:n);
    H(k+1:n, k:n) = H(k+1:n, k:n) - w * v;
    v = H(1:n, k+1:n) * w;
    H(1:n, k+1:n) = H(1:n, k+1:n) - v * w';
    H(k+2:n, k) = 0;
    if nargout == 2
        U(k+1:n,k) = w;
    end
end

if nargout == 2
    U(:, n) = zeros(n, 1); U(n, n) = 1;
    U(:, n-1) = zeros(n, 1); U(n-1, n-1) = 1;
    for k = n-2:-1:1
        w = U(k+1:n, k);
        v = w' * U(k+1:n, k+1:n);
        U(k+1:n, k+1:n) = U(k+1:n, k+1:n) - w * v;
        U(:, k) = zeros(n, 1); U(k, k) = 1;
    end
end