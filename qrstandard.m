function [H, Q] = qrstandard(A, shift, tol)
% QRSTANDARD    Compute Schur decomposition of a square matrix
%
% Given a square matrix A, this function computes its Schur decomposition using
% single-shifted QR algorithm.
%
% argin:
%   A     - A square matrix
%   shift - A string of shift type (options: 'none', 'rayleigh' or 'wilkinson';
%           default:'wilkinson')
%   tol   - Tolerance of the precision (default: eps)
%
% argout:
%   T - The upper triangular matrix of Schur decomposition Q'AQ = T
%   Q - The unitray matrix of Schur decomposition Q'AQ = T
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-03-31
% -------------------------------------------------

if nargin == 1
    shift = 'wilkinson';
    tol = eps;
elseif nargin == 2
    tol = eps;
end

n = size(A);
if n(1) ~= n(2)
    error('Fatal Error: The input is not a sqaure matrix!');
end
n = n(1);
if n == 1
    H = A;
    Q = 1;
    return
end

% Step 1: Upper Hessenberg reduction
[H, U] = qrstandard_reduction(A);
if nargout == 2
    Q = U;
end

m = 0;
while true
    % Step 2: Judge convergence of an eigenvalue
    for k = 1:n-1
        if abs(H(k+1, k)) <= tol * (abs(H(k, k)) + abs(H(k+1, k+1)))
            H(k+1, k) = 0;
        end
    end
    for m = m:n-2
        if H(n-m, n-m-1) ~=0
            break
        end
    end
    if m == n-2 && H(n-m, n-m-1) == 0  % Final condition
        break
    end
    l = 0;
    for k = n-m-2:-1:1
        if H(k+1, k) == 0
            l = k;
            break
        end
    end
    
    % Step 3: QR iteration
    [H22, G] = qrstandard_iteration(H(l+1:n-m, l+1:n-m), shift);
    H(l+1:n-m, l+1:n-m) = H22;
    for k = l+1:n-m-1
        H(1:l, k:k+1) = H(1:l, k:k+1) * G(:, :, k-l)';
        H(k:k+1, n-m+1:n) = G(:, :, k-l) * H(k:k+1, n-m+1:n);
        if nargout == 2
            Q(:, k:k+1) = Q(:, k:k+1) * G(:, :, k-l)';
        end
    end
end