function [E, A, Q, Z] = qzperiodic_iteration(E, A, l, m, shift, Q, Z)
% QZPERIODIC_ITERATION    Perform a single step of QZ iteration for 
%                         periodic eigenvalue problem.
%
% This version of QZ iteration performs on periodic coefficients (E{j}, A{j}),
% j = 1:n, in Hessenberg-triangular form, using double-shift strategy if needed.
%
% argin:
%   E, A  - Two cell arrays of coefficients, both with the length of p, where
%           all E matrices are upper triangular, all A matrices but A{1} are
%           upper triangular, and A{1} is upper-Hessenberg.
%   l, m  - Two block-deflation value in QZ iteration.
%   shift - Shift type ('double' or 'comb', default: 'comb').
%   Q, Z  - Current cell arrays of unitray matrices Q and Z to be accumulated in
%           this iteration step.
%
% argout:
%   E, A, Q, Z - The cell arrays of matrices where a QZ iteration step have been
%                performed. Note that Q and Z would be computed only when
%                needed.
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-04-22
% -------------------------------------------------

% Check inputs
p = length(E);
n = length(E{1});
n2 = n - m;

if n2 == 1
    return;
end

% Compute the Wilkinson's shift
H = eye(3, 3);
for j = p:-1:2
    H = H / E{j}(n2-2:n2, n2-2:n2);
    H = H * A{j}(n2-2:n2, n2-2:n2);
end
H = A{1}(n2-2:n2, n2-2:n2) * H / E{1}(n2-2:n2, n2-2:n2);
mu = eig(H(2:3, 2:3));

% QZ iteration
if isreal(mu) && isequal(shift, 'comb')  % Single-shift implicit
    
    % Implicit Q: The first column
    if abs(mu(1) - H(3, 3)) < abs(mu(2) - H(3, 3))
        mu = mu(1);
    else
        mu = mu(2);
    end
    mu = mu * E{1}(l+1, l+1);
    for j = 2:p
        mu = mu * E{j}(l+1, l+1) / A{j}(l+1, l+1);
    end
    a = A{1}(l+1:l+2, l+1) - [mu; 0];
    G = givens(a(1), a(2));
    A{1}(l+1:l+2, l+1:n) = G * A{1}(l+1:l+2, l+1:n);
    E{1}(l+1:l+2, l+1:n) = G * E{1}(l+1:l+2, l+1:n);
    if nargout >= 3
        Q{1}(1:n, l+1:l+2) = Q{1}(1:n, l+1:l+2) * G';
    end
    
    % Iteration
    for i = l+1:n2-1
        for j = 1:p-1
            % On E{1:p-1}
            G = givens(-E{j}(i+1, i+1), E{j}(i+1, i));
            E{j}(1:i+1, i:i+1) = E{j}(1:i+1, i:i+1) * G';
            A{j+1}(1:i+1, i:i+1) =  A{j+1}(1:i+1, i:i+1) * G';
            if nargout == 4
                Z{j}(1:n, i:i+1) = Z{j}(1:n, i:i+1) * G';
            end
            
            % On A{2:p}
            G = givens(A{j+1}(i, i), A{j+1}(i+1, i));
            A{j+1}(i:i+1, i:n) = G * A{j+1}(i:i+1, i:n);
            E{j+1}(i:i+1, i:n) = G * E{j+1}(i:i+1, i:n);
            if nargout >= 3
                Q{j+1}(1:n, i:i+1) = Q{j+1}(1:n, i:i+1) * G';
            end
        end
        
        % On E{p} and A{1}
        G = givens(-E{p}(i+1, i+1), E{p}(i+1, i));
        E{p}(1:i+1, i:i+1) = E{p}(1:i+1, i:i+1) * G';
        A{1}(1:n2, i:i+1) =  A{1}(1:n2, i:i+1) * G';
        if nargout == 4
            Z{p}(1:n, i:i+1) = Z{p}(1:n, i:i+1) * G';
        end
        if i < n2-1
            G = givens(A{1}(i+1, i), A{1}(i+2, i));
            A{1}(i+1:i+2, i:n) = G * A{1}(i+1:i+2, i:n);
            E{1}(i+1:i+2, i:n) = G * E{1}(i+1:i+2, i:n);
            if nargout >= 3
                Q{1}(1:n, i+1:i+2) = Q{1}(1:n, i+1:i+2) * G';
            end
        end
    end

else  % Double-shift implicit
    
    % Implicit Q: The first column
    s = sum(mu);
    t = prod(mu);
    
    H = A{1}(l+1:l+3, l+1:l+2);
    beta = 1;
    for j = p:-1:1
        H = H / E{j}(l+1:l+2, l+1:l+2);
        H = H * A{j}(l+1:l+2, l+1:l+2);
        if j == 1
            beta = beta * E{j}(l+1, l+1);
        else
            beta = beta * E{j}(l+1, l+1) / A{j}(l+1, l+1);
        end
    end
    a = H(1:3, 1) - s*A{1}(l+1:l+3, l+1) + [t; 0; 0]*beta;
    
    w = house(a);
    v = w' * A{1}(l+1:l+3, l+1:n);
    A{1}(l+1:l+3, l+1:n) = A{1}(l+1:l+3, l+1:n) - w * v;
    v = w' * E{1}(l+1:l+3, l+1:n);
    E{1}(l+1:l+3, l+1:n) = E{1}(l+1:l+3, l+1:n) - w * v;
    if nargout >= 3
        v = Q{1}(1:n, l+1:l+3) * w;
        Q{1}(1:n, l+1:l+3) = Q{1}(1:n, l+1:l+3) - v * w';
    end
    
    % Iteration
    for i = l+1:n2-2
        for j = 1:p-1
            % On E{1:p-1}
            w = house(E{j}(i+2, i:i+2));
            H = flip(eye(3)) - w * flip(w');
            E{j}(1:i+2, i:i+2) = E{j}(1:i+2, i:i+2) * H;
            A{j+1}(1:i+2, i:i+2) = A{j+1}(1:i+2, i:i+2) * H;
            G = givens(-E{j}(i+1, i+1), E{j}(i+1, i));
            E{j}(1:i+1, i:i+1) = E{j}(1:i+1, i:i+1) * G';
            A{j+1}(1:i+2, i:i+1) = A{j+1}(1:i+2, i:i+1) * G';
            if nargout == 4
                Z{j}(1:n, i:i+2) = Z{j}(1:n, i:i+2) * H;
                Z{j}(1:n, i:i+1) = Z{j}(1:n, i:i+1) * G';
            end
            
            % On A{2:p}
            w = house(A{j+1}(i:i+2, i));
            v = w' * A{j+1}(i:i+2, i:n);
            A{j+1}(i:i+2, i:n) = A{j+1}(i:i+2, i:n) - w * v;
            v = w' * E{j+1}(i:i+2, i:n);
            E{j+1}(i:i+2, i:n) = E{j+1}(i:i+2, i:n) - w * v;
            G = givens(A{j+1}(i+1, i+1), A{j+1}(i+2, i+1));
            A{j+1}(i+1:i+2, i:n) = G * A{j+1}(i+1:i+2, i:n);
            E{j+1}(i+1:i+2, i:n) = G * E{j+1}(i+1:i+2, i:n);
            if nargout >= 3
                v = Q{j+1}(1:n, i:i+2) * w;
                Q{j+1}(1:n, i:i+2) = Q{j+1}(1:n, i:i+2) - v * w';
                Q{j+1}(1:n, i+1:i+2) = Q{j+1}(1:n, i+1:i+2) * G';
            end
        end
        
        % On E{p} and A{1}
        w = house(E{p}(i+2, i:i+2));
        H = flip(eye(3)) - w * flip(w');
        E{p}(1:i+2, i:i+2) = E{p}(1:i+2, i:i+2) * H;
        A{1}(1:n2, i:i+2) = A{1}(1:n2, i:i+2) * H;
        G = givens(-E{p}(i+1, i+1), E{p}(i+1, i));
        E{p}(1:i+1, i:i+1) = E{p}(1:i+1, i:i+1) * G';
        A{1}(1:n2, i:i+1) = A{1}(1:n2, i:i+1) * G';
        if nargout == 4
            Z{p}(1:n, i:i+2) = Z{p}(1:n, i:i+2) * H;
            Z{p}(1:n, i:i+1) = Z{p}(1:n, i:i+1) * G';
        end
        
        if i < n2-2
            w = house(A{1}(i+1:i+3, i));
            v = w' * A{1}(i+1:i+3, i:n);
            A{1}(i+1:i+3, i:n) = A{1}(i+1:i+3, i:n) - w * v;
            v = w' * E{1}(i+1:i+3, i:n);
            E{1}(i+1:i+3, i:n) = E{1}(i+1:i+3, i:n) - w * v;
            if nargout >= 3
                v = Q{1}(1:n, i+1:i+3) * w;
                Q{1}(1:n, i+1:i+3) = Q{1}(1:n, i+1:i+3) - v * w';
            end
        end
    end
    
    % Final step
    for j = 1:p
        if j == 1
            G = givens(A{j}(n2-1, n2-2), A{j}(n2, n2-2));
            A{j}(n2-1:n2, n2-2:n) = G * A{j}(n2-1:n2, n2-2:n);
            E{j}(n2-1:n2, n2-1:n) = G * E{j}(n2-1:n2, n2-1:n);
            if nargout >= 3
                Q{j}(1:n, n2-1:n2) = Q{j}(1:n, n2-1:n2) * G';
            end
        else
            G = givens(A{j}(n2-1, n2-1), A{j}(n2, n2-1));
            A{j}(n2-1:n2, n2-1:n) = G * A{j}(n2-1:n2, n2-1:n);
            E{j}(n2-1:n2, n2-1:n) = G * E{j}(n2-1:n2, n2-1:n);
            if nargout >= 3
                Q{j}(1:n, n2-1:n2) = Q{j}(1:n, n2-1:n2) * G';
            end
        end
        
        G = givens(-E{j}(n2, n2), E{j}(n2, n2-1));
        E{j}(1:n2, n2-1:n2) = E{j}(1:n2, n2-1:n2) * G';
        if j < p
            A{j+1}(1:n2, n2-1:n2) = A{j+1}(1:n2, n2-1:n2) * G';
        else
            A{1}(1:n2, n2-1:n2) = A{1}(1:n2, n2-1:n2) * G';
        end
        if nargout == 4
            Z{j}(1:n, n2-1:n2) = Z{j}(1:n, n2-1:n2) * G';
        end
    end
    
end

% Set zeros
for j = 1:p
    if j == 1
        A{j} = triu(A{j}, -1);
    else
        A{j} = triu(A{j});
    end
    E{j} = triu(E{j});
end