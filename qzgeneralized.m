function [A, B, Q, Z] = qzgeneralized(A, B, shift, tol)
% QZGENERALIZED    Compute generalized Schur decomposition of a matrix pair
%                  (A, B).
%
% Given a real matrix pair (A, B), this function computes its generalized Schur
% decomposition using QZ algorithm.
%
% argin:
%   A, B  - The real matrix pair (A, B) to compute generalized Schur
%           decomposition.
%   shift - String of shift type ('single', 'double' or 'comb', default: 'comb).
%   tol   - Tolerance of the precision (default: eps).
%
% argout:
%   S, T, Q, Z - Upper triangular matrices and unitray matrices of the
%                generalized Schur decomposition, such that:
%                Q'AZ = S and Q'BZ = T.
%                Note that Q and Z would be computed only when needed.
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-04-14
% -------------------------------------------------

if nargin < 3
    shift = 'comb';
    tol = eps;
elseif nargin < 4
    tol = eps;
end

n = size(A); nb = size(B);
if n(1) ~= n(2) || n(2) ~= nb(1) || nb(1) ~= nb(2)
    error('Fatal Error: The input is not an appropriate matrix pair!');
end
n = n(1);

if n == 1
    Q = 1;
    Z = 1;
    return
end

% Step 1: Hessenberg-triangular reduction
if nargout <= 2
    [A, B] = qzgeneralized_reduction(A, B);
elseif nargout == 3
    [A, B, Q] = qzgeneralized_reduction(A, B);
elseif nargout == 4
    [A, B, Q, Z] = qzgeneralized_reduction(A, B);
end

if n == 2
    return
end

m = 0;
while true
    % Step 2: Convergence judgement and block deflation
    for k = 1:n-1
        if abs(A(k+1, k)) <= tol * (abs(A(k, k)) + abs(A(k+1, k+1)))
            A(k+1, k) = 0;
        end
    end
    for m = m:n-3
        if A(n-m, n-m-1) ~= 0 && A(n-m-1, n-m-2) ~= 0
            break
        end
    end
    if (m == n-3) && (A(2, 1) == 0 || A(3, 2) == 0)  % Final condition
        break
    end
    l = 0;
    for k = n-m-2:-1:1
        if A(k+1,k) == 0
            l = k;
            break
        end
    end
    % Step 3: QZ deflation
    if nargout <= 2
        [A, B] = qzgeneralized_deflation(A, B, l, m);
    elseif nargout == 3
        [A, B, Q] = qzgeneralized_deflation(A, B, l, m, Q);
    elseif nargout == 4
        [A, B, Q, Z] = qzgeneralized_deflation(A, B, l, m, Q, Z);
    end
    for m = m:n-3
        if A(n-m, n-m-1) ~= 0 && A(n-m-1, n-m-2) ~= 0
            break
        end
    end
    
    % Step 4: QZ iteration
    if nargout <= 2
        [A, B] = qzgeneralized_iteration(A, B, l, m, shift);
    elseif nargout == 3
        [A, B, Q] = qzgeneralized_iteration(A, B, l, m, shift, Q);
    elseif nargout == 4
        [A, B, Q, Z] = qzgeneralized_iteration(A, B, l, m, shift, Q, Z);
    end
end