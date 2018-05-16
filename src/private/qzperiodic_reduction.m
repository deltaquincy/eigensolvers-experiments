function [E, A, Q, Z] = qzperiodic_reduction(E, A)
% QZPERIODIC_REDUCTION    Compute Hessenberg-triangular form of coefficients of
%                         a periodic system (E{j}, A{j}), j = 1:p.
%
% Given coefficients of a periodic system (E{j}, A{j}), j = 1:p, this function 
% computes the upper Hessenberg-triangular form of it, where A{1} is upper
% Hessenberg, and the others are upper triangular.
%
% argin:
%   E, A - Two cell arrays of coefficients, both with the length of p.
%
% argout:
%   E, A - Transformed cell arrays of coefficients, where A{1} is upper
%          Hessenberg, and the others are upper triangular.
%   Q, Z - Cell arrays of unitray matrices, so that, Q{j}'*E{j}*Z{j} are upper
%          triangular, Q{j}'*A{j}*Z{j-1} are also upper triangular but j = 1.
%          When j = 1, Q{1}'*A{1}*Z{p} is upper Hessenberg.
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-04-21
% -------------------------------------------------

% Check inputs
p = length(E);
n = length(E{1});

% Initialize unitray matrices
if nargout == 3
    Q = cell(1, p);
elseif nargout == 4
    Q = cell(1, p);
    Z = cell(1, p);
end

% QR decomposition of E matrices
if nargout < 3
    for j = 1:p
        [U, E{j}] = qr(E{j});
        A{j} = U * A{j};
    end
else
    for j = 1:p
        Z{j} = eye(n, n);
        [Q{j}, E{j}] = qr(E{j});
        A{j} = Q{j}' * A{j};
    end
end

% Reduction of A matrices
for i = 1:n-1
    % Upper triangular for A{2:end}
    for j = 2:p
        for k = n:-1:i+1                
            G = givens(A{j}(k-1, i), A{j}(k, i));
            A{j}(k-1:k, i:n) = G * A{j}(k-1:k, i:n);
            E{j}(k-1:k, k-1:n) = G * E{j}(k-1:k, k-1:n);
            if nargout >= 3
                Q{j}(1:n, k-1:k) = Q{j}(1:n, k-1:k) * G';
            end
            
            G = givens(-E{j}(k, k), E{j}(k, k-1));
            E{j}(1:k, k-1:k) = E{j}(1:k, k-1:k) * G';
            if j < p
                A{j+1}(1:n, k-1:k) = A{j+1}(1:n, k-1:k) * G';
            else
                A{1}(1:n, k-1:k) = A{1}(1:n, k-1:k) * G';
            end
            if nargout == 4
                Z{j}(1:n, k-1:k) = Z{j}(1:n, k-1:k) * G';
            end
        end
    end
    
    % Upper Hessenberg for A{1}
    for k = n:-1:i+2
        G = givens(A{1}(k-1, i), A{1}(k, i));
        A{1}(k-1:k, i:n) = G * A{1}(k-1:k, i:n);
        E{1}(k-1:k, k-1:n) = G * E{1}(k-1:k, k-1:n);
        if nargout >= 3
            Q{1}(1:n, k-1:k) = Q{1}(1:n, k-1:k) * G';
        end
        
        G = givens(-E{1}(k, k), E{1}(k, k-1));
        E{1}(1:k, k-1:k) = E{1}(1:k, k-1:k) * G';
        A{2}(1:n, k-1:k) = A{2}(1:n, k-1:k) * G';
        if nargout == 4
            Z{1}(1:n, k-1:k) = Z{1}(1:n, k-1:k) * G';
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
