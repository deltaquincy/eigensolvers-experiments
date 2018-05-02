function [A, B, Q, Z] = qzgeneralized_iteration(A, B, l, m, shift, Q, Z)
% QZGENERALIZED_ITERATION    Perform a single step of QZ iteration for 
%                            generalized eigenvalue problem.
%
% This version of QZ iteration performs on a matrix pair (A, B) in 
% Hessenberg-triangular form, using double-shift strategy if needed.
%
% argin:
%   A, B  - A matrix pair (A, B) in Hessenberg-triangular form to perform QZ
%           iteration step.
%   l, m  - Two block-deflation value in QZ iteration.
%   shift - String of shift type ('single', 'double' or 'comb', default: 'comb).
%   Q, Z  - Current unitray matrices Q and Z to be accumulated in this iteration
%           step.
%
% argout:
%   A, B, Q, Z - The matrices where a QZ iteration step have been performed.
%                Note that Q and Z would be computed only when needed.
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-04-14
% -------------------------------------------------

n = length(A);
n2 = n-m;

mu = eig(A(n2-1:n2, n2-1:n2), B(n2-1:n2, n2-1:n2));

if isequal(shift, 'single') || (isequal(shift, 'comb') && isreal(mu))
    % Single-shift QZ iteration
    
    % TODO: Wilkinson shift
    mudist = abs(mu - A(n2, n2)/B(n2, n2));
    if mudist(1) <= mudist(2)
        mu = mu(1);
    else
        mu = mu(2);
    end
    
    GZ = zeros(2, 2, n2-1);
    A = A - mu * B;
    for k = l+1:n2-1
        [c, s] = givens(A(k, k), A(k+1, k));
        G = [c, s; -conj(s), conj(c)];
        A(k:k+1, k:n) = G * A(k:k+1, k:n);
        B(k:k+1, k:n) = G * B(k:k+1, k:n);
        if nargout >= 3
            Q(1:n, k:k+1) = Q(1:n, k:k+1) * G';
        end
        [c, s] = givens(-B(k+1, k+1), B(k+1, k));
        GZ(:, :, k-l) = [c, s; -conj(s), conj(c)];
        B(1:k+1, k:k+1) = B(1:k+1, k:k+1) * GZ(:, :, k-l)';
        B(k+1, k) = 0;
        if nargout == 4
            Z(1:n, k:k+1) = Z(1:n, k:k+1) * GZ(:, :, k-l)';
        end
    end
    for k = l+1:n2-1
        A(1:k+1, k:k+1) = A(1:k+1, k:k+1) * GZ(:, :, k-l)';
    end
    A = A + mu * B;
    
else
    % Double-shift implicit QZ iteration
    s = sum(mu);
    t = prod(mu);
    H = A(l+1:l+3, l+1:l+3) / B(l+1:l+3, l+1:l+3);
    
    xi1 = H(1, 1)^2 + H(1, 2)*H(2, 1) - s*H(1, 1) + t;
    xi2 = H(2, 1) * (H(1, 1) + H(2, 2) - s);
    xi3 = H(2, 1) * H(3, 2);
    

    for k = l+1:n2-2
        % Q step
        w = house([xi1; xi2; xi3]);
        v = w' * A(k:k+2, l+1:n);  % TODO
        A(k:k+2, l+1:n) = A(k:k+2, l+1:n) - w * v;
        v = w' * B(k:k+2, k:n);
        B(k:k+2, k:n) = B(k:k+2, k:n) - w * v;
        if nargout >= 3
            v = Q(1:n, k:k+2) * w;
            Q(1:n, k:k+2) = Q(1:n, k:k+2) - v * w';
        end

        % Z-step
        w = house(B(k+2, k:k+2));  % TODO
        H = flip(eye(3)) - w * flip(w');
        A(1:n2, k:k+2) = A(1:n2, k:k+2) * H;
        B(1:k+2, k:k+2) = B(1:k+2, k:k+2) * H;
        [c, s] = givens(-B(k+1, k+1), B(k+1, k));
        G = [c, s; -conj(s), conj(c)];
        A(1:n2, k:k+1) = A(1:n2, k:k+1) * G';
        B(1:k+1, k:k+1) = B(1:k+1, k:k+1) * G';
        if nargout == 4
            Z(1:n, k:k+2) = Z(1:n, k:k+2) * H;
            Z(1:n, k:k+1) = Z(1:n, k:k+1) * G';
        end

        % Update xi1, xi2 and xi3
        xi1 = A(k+1, k);
        xi2 = A(k+2, k);
        if k < n-2
            xi3 = A(k+3, k);
        end
    end

    [c, s] = givens(xi1, xi2);
    G = [c, s; -conj(s), conj(c)];
    A(n2-1:n2, n2-2:n) = G * A(n2-1:n2, n2-2:n);
    B(n2-1:n2, n2-1:n) = G * B(n2-1:n2, n2-1:n);
    if nargout >= 3
        Q(1:n, n2-1:n2) = Q(1:n, n2-1:n2) * G';
    end

    [c, s] = givens(-B(n2, n2), B(n2, n2-1));
    G = [c, s; -conj(s), conj(c)];
    A(1:n2, n2-1:n2) = A(1:n2, n2-1:n2) * G';
    B(1:n2, n2-1:n2) = B(1:n2, n2-1:n2) * G';
    if nargout == 4
        Z(1:n, n2-1:n2) = Z(1:n, n2-1:n2) * G';
    end

    A = triu(A, -1);
    B = triu(B);
    
end