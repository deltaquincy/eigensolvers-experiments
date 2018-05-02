function [E, A, Q, Z] = qzperiodic_deflation(E, A, l, m, Q, Z)
% QZPERIODIC_DEFLATION    Perform deflation strategy in periodic QZ algorithm.
%
% Given Reducted coefficients of a periodic system (E{j}, A{j}), j = 1:p, this
% function performs deflation strategy when any of E matrices or A matrices is
% singular.
%
% argin:
%   E, A - Two cell arrays of coefficients, both with the length of p, where all
%          E matrices are upper triangular, all A matrices but A{1} are
%          upper triangular, and A{1} is upper-Hessenberg.
%   l, m - Two block-deflation value in QZ iteration.
%   Q, Z - Current unitray matrices Q and Z to be accumulated in this iteration
%          step.
%
% argout:
%   E, A - Transformed cell arrays of coefficients.
%   Q, Z - Cell arrays of unitray matrices.
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

tol = 1e-7;


% -----------------------
% Deflation on E matrices
% -----------------------

jj = 0;
ii = 0;

% Check diagonal elements
for j = 1:p
    for i = l+1:n2
        if abs(E{j}(i, i)) < tol
            jj = j;
            ii = i;
            break;
        end
    end
end

if jj == 0 && ii <= 2  % TODO
    return
end

for i = ii:n-1
    for j = jj:-1:1
        if j == jj
            G = givens(E{j}(i, i+1), E{j}(i+1, i+1));
        else
            G = givens(E{j}(i, i), E{j}(i+1, i));
        end
        E{j}(i:i+1, i:n) = G * E{j}(i:i+1, i:n);
        if j == 1
            A{1}(i:i+1, i-1:n) = G * A{1}(i:i+1, i-1:n);
            if nargout >= 3
                Q{j}(1:n, i:i+1) = Q{j}(1:n, i:i+1) * G';
            end
            G = givens(-A{1}(i+1, i), A{1}(i+1, i-1));
            A{1}(1:n, i-1:i) = A{1}(1:n, i-1:i) * G';
            E{p}(1:n, i-1:i) = E{p}(1:n, i-1:i) * G';
            if nargout == 4
                Z{p}(1:n, i-1:i) = Z{p}(1:n, i-1:i) * G';
            end
        else
            A{j}(i:i+1, i:n) = G * A{j}(i:i+1, i:n);
            if nargout >= 3
                Q{j}(1:n, i:i+1) = Q{j}(1:n, i:i+1) * G';
            end
            G = givens(-A{j}(i+1, i+1), A{j}(i+1, i));
            A{j}(1:n, i:i+1) = A{j}(1:n, i:i+1) * G';
            E{j-1}(1:n, i:i+1) = E{j-1}(1:n, i:i+1) * G';
            if nargout == 4
                Z{j-1}(1:n, i:i+1) = Z{j-1}(1:n, i:i+1) * G';
            end
        end
    end
    
    for j = p:-1:jj+1
        G = givens(E{j}(i-1, i-1), E{j}(i, i-1));
        E{j}(i-1:i, 1:n) = G * E{j}(i-1:i, 1:n);
        A{j}(i-1:i, 1:n) = G * A{j}(i-1:i, 1:n);
        if nargout >= 3
            Q{j}(1:n, i-1:i) = Q{j}(1:n, i-1:i) * G';
        end
        
        G = givens(-A{j}(i, i), A{j}(i, i-1));
        A{j}(1:n, i-1:i) = A{j}(1:n, i-1:i) * G';
        E{j-1}(1:n, i-1:i) = E{j-1}(1:n, i-1:i) * G';
        if nargout == 4
            Z{j-1}(1:n, i-1:i) = Z{j-1}(1:n, i-1:i) * G';
        end
    end
end

G = givens(-A{1}(n, n), A{1}(n, n-1));
A{1}(1:n, n-1:n) = A{1}(1:n, n-1:n) * G';
E{p}(1:n, n-1:n) = E{p}(1:n, n-1:n) * G';
if nargout == 4
    Z{p}(1:n, n-1:n) = Z{p}(1:n, n-1:n) * G';
end
for j = p:-1:jj+1
    G = givens(E{j}(n-1, n-1), E{j}(n, n-1));
    E{j}(n-1:n, 1:n) = G * E{j}(n-1:n, 1:n);
    A{j}(n-1:n, 1:n) = G * A{j}(n-1:n, 1:n);
    if nargout >= 3
        Q{j}(1:n, n-1:n) = Q{j}(1:n, n-1:n) * G';
    end
    
    G = givens(-A{j}(n, n), A{j}(n, n-1));
    A{j}(1:n, n-1:n) = A{j}(1:n, n-1:n) * G';
    E{j-1}(1:n, n-1:n) = E{j-1}(1:n, n-1:n) * G';
    if nargout == 4
        Z{j-1}(1:n, n-1:n) = Z{j-1}(1:n, n-1:n) * G';
    end
end


% -----------------------
% Deflation on A matrices
% -----------------------

% *** Unfinished ***

% jj = 0;
% ii = 0;
% 
% % Check diagonal elements
% for j = 2:p
%     for i = l+1:n2
%         if abs(A{j}(i, i)) < tol
%             jj = j;
%             ii = i;
%             break;
%         end
%     end
% end
% 
% if jj == 0 && ii <= 1  % TODO
%     return
% end

% for i = 1:ii-1
%     for j = 1:jj
%         G = givens(A{j}(i, i), A{j}(i+1, i));
%         A{j}(i:i+1, 1:n) = G * A{j}(i+i+1, 1:n);
%         E{j}(i:i+1, 1:n) = G * E{j}(i:i+1, 1:n);
%         if nargout >= 3
%             Q{j}(1:n, i:i+1) = Q{j}(1:n, i:i+1) * G';
%         end
%         
%         G = givens(-E{j}(i+1, i+1), E{j}(i+1, i));
%         E{j}(1:n, i:i+1) = E{j}(1:n, i:i+1) * G';
%         A{j+1}(1:n, i:i+1) = A{j+1}(1:n, i:i+1) * G';
%         if nargout == 4
%             Z{j}(1:n, i:i+1) = Z{j}(1:n, i:i+1) * G';
%         end
%     end
%     
% end
    
