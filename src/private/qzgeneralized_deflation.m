function [A, B, Q, Z] = qzgeneralized_deflation(A, B, l, m, Q, Z)
% QZGENERALIZED_DEFLATION    Perform deflation strategy in QZ algorithm.
%
% In the process of QZ algorithm, this function performs the third step:
% deflation on the singular matrix B.
%
% argin:
%   A, B - A matrix pair (A, B) in Hessenberg-triangular form to perform QZ
%          deflation step.
%   l, m - Two block-deflation value in QZ iteration.
%   Q, Z - Current unitray matrices Q and Z to be accumulated in this iteration
%          step.
%
% argout:
%   A, B, Q, Z - The matrices where QZ iteration step have been performed QZ
%                deflation.
%                Note that Q and Z would be computed only when needed.
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-04-14
% -------------------------------------------------

n = length(A);
tol = 1e-10;

z = 0;  % The number of zero elements
while true
    % Check whether zero element exists
    for k = l+1:n-m-z
        % TODO: tolerance?
        if k == n
            tol = abs(B(k-1, k)) * tol;
        else
            tol = abs(B(k, k+1)) * tol;
        end
        if abs(B(k, k)) < tol
            B(k, k) = 0;
            z = z + 1;
            break
        end
    end
    if k == n-m-z  % Final condition
        return
    end
    if k == 1
        [c, s] = givens(B(1, 2), B(2, 2));
        G = [c, s; -conj(s), conj(c)];
        A(1:2, 1:n) = G * A(1:2, 1:n);
        B(1:2, 1:n) = G * B(1:2, 1:n);
        B(2, 2) = 0;
        if nargout >= 3
            Q(1:n, 1:2) = Q(1:n, 1:2) * G';
        end
        [c, s] = givens(B(2, 3), B(3, 3));
        G = [c, s; -conj(s), conj(c)];
        A(2:3, 1:n) = G * A(2:3, 1:n);
        B(2:3, 1:n) = G * B(2:3, 1:n);
        B(3, 3) = 0;
        if nargout >= 3
            Q(1:n, 2:3) = Q(1:n, 2:3) * G';
        end
        [c, s] = givens(-A(3, 2), A(3, 1));
        G = [c, s; -conj(s), conj(c)];
        A(1:3, 1:2) = A(1:3, 1:2) * G';
        B(1:3, 1:2) = B(1:3, 1:2) * G';
        k = k + 1;
    end
    
    for i = k:n-1
        [c, s] = givens(B(i, i+1), B(i+1, i+1));
        G = [c, s; -conj(s), conj(c)];
        A(i:i+1, i-1:n) = G * A(i:i+1, i-1:n);
        B(i:i+1, i+1:n) = G * B(i:i+1, i+1:n);
        B(i+1, i+1) = 0;
        if nargout >= 3
            Q(1:n, i:i+1) = Q(1:n, i:i+1) * G';
        end

        [c, s] = givens(-A(i+1, i), A(i+1, i-1));
        G = [c, s; -conj(s), conj(c)];
        A(1:i+1, i-1:i) = A(1:i+1, i-1:i) * G';
        A(i+1, i-1) = 0;
        B(1:i-1, i-1:i) = B(1:i-1, i-1:i) * G';
        if nargout == 4
            Z(1:n, i-1:i) = Z(1:n, i-1:i) * G';
        end
    end

    [c, s] = givens(-A(n, n), A(n, n-1));
    G = [c, s; -conj(s), conj(c)];
    A(1:n, n-1:n) = A(1:n, n-1:n) * G';
    A(n, n-1) = 0;
    B(1:n-1, n-1:n) = B(1:n-1, n-1:n) * G';
    if nargout == 4
        Z(1:n, n-1:n) = Z(1:n, n-1:n) * G';
    end
end



