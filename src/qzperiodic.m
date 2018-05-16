function [E, A, Q, Z, iter] = qzperiodic(E, A, shift, tol)
% QZPERIODIC    Compute periodic Schur decomposition of coefficients of a linear
%               system (E{j}, A{j}) where j = 1:n.
%
% Given coefficients of a periodic system (E{j}, A{j}), j = 1:p, this function 
% computes Schur canonical form of the periodic system, where all E and A
% matrices but A{1} are upper triangular, and A{1} is quasi-triangular.
%
% argin:
%   E, A  - Two cell arrays of coefficients, both with the length of p.
%   shift - Shift type ('double' or 'comb', default: 'comb').
%   tol   - Tolerance of the precision (default: 1e-10).
%
% argout:
%   E, A - Transformed cell arrays of coefficients, where A{1} is
%          quasi-triangular, and the others are upper triangular.
%   Q, Z - Cell arrays of unitray matrices, so that, Q{j}'*E{j}*Z{j} are upper
%          triangular, Q{j}'*A{j}*Z{j-1} are also upper triangular but j = 1.
%          When j = 1, Q{1}'*A{1}*Z{p} is quasi-triangular.
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-04-22
% -------------------------------------------------

if nargin < 3
    shift = 'comb';
    tol = 1e-10;
elseif nargin < 4
    tol = 1e-10;
end

% Check input
p = length(E);
if length(A) ~= p
    error('Fatal Error: The input is not approriate.');
end
n = size(E{1});
for j = 1:p
    if ~isequal(size(E{j}), n) || ~isequal(size(A{j}), n)
        error('Fatal Error: The input is not approriate.');
    end
end
n = n(1);

% Step 1: Reduction
if nargout <= 2
    [E, A] = qzperiodic_reduction(E, A);
elseif nargout == 3
    [E, A, Q] = qzperiodic_reduction(E, A);
elseif nargout >= 4
    [E, A, Q, Z] = qzperiodic_reduction(E, A);
end

iter = 0;
while true
    % Step 2: Judge convergence
    for i = 1:n-1
        if abs(A{1}(i+1, i)) <= tol * (abs(A{1}(i, i)) + abs(A{1}(i+1, i+1)))
            A{1}(i+1, i) = 0;
        end
    end
    [l, m] = qzperiodic_judge(E, A);
    if m == n
        break;
    end
    
    % Step 3: Deflation and judge convergence again
    if nargout <= 2
        [E, A] = qzperiodic_deflation(E, A, l, m);
    elseif nargout == 3
        [E, A, Q] = qzperiodic_deflation(E, A, l, m, Q);
    elseif nargout >= 4
        [E, A, Q, Z] = qzperiodic_deflation(E, A, l, m, Q, Z);
    end
    [l, m] = qzperiodic_judge(E, A);
    if m == n
        break;
    end
    
    % Step 4: QZ iteration
    if nargout <= 2
        [E, A] = qzperiodic_iteration(E, A, l, m, shift);
    elseif nargout == 3
        [E, A, Q] = qzperiodic_iteration(E, A, l, m, shift, Q);
    elseif nargout >= 4
        [E, A, Q, Z] = qzperiodic_iteration(E, A, l, m, shift, Q, Z);
    end
    iter = iter + 1;
end