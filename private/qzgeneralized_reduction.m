function [A, B, Q, Z] = qzgeneralized_reduction(A, B, anm)
% QZGENERALIZED_REDUCTION    Compute Hessenberg-triangular form of a matrix
%                            pair (A, B).
%
% Given a matrix pair (A, B), this function computes an upper Hessenberg matrix
% Q'AZ and an upper triangular matrix Q'BZ, where Q and Z are unitray.
%
% argin:
%   A, B - The matrix pair to perform Hessenberg-triangular reduction.
%
% argout:
%   A, B, Q, Z - The matrices in Hessenberg-triangular form, such that A is an
%                upper Hessenberg form of the original A, B an upper triangular
%                form of the original B, and Q, Z are unitray matrices such that
%                A = Q'AZ and B = Q'BZ.
%                Note that Q and Z would be computed only when needed.
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-04-11
% -------------------------------------------------

if nargin < 3
    anm = false;
end    

n = size(A); nb = size(B);
if n(1) ~= n(2) || n(2) ~= nb(1) || nb(1) ~= nb(2)
    error('Fatal Error: The input is not an appropriate matrix pair!');
end
n = n(1);

if anm
    subplot(1, 2, 1);
    heatmap(abs(A)); title('A'); getframe(gcf);
    subplot(1, 2, 2);
    heatmap(abs(B)); title('B'); getframe(gcf);
end

[Q, B] = qr(B);  % Here we use the built-in QR decomposition function
A = Q' * A;

if anm
    subplot(1, 2, 1);
    heatmap(abs(A)); title('A'); getframe(gcf);
    subplot(1, 2, 2);
    heatmap(abs(B)); title('B'); getframe(gcf);
end

if nargout == 4
    Z = eye(n, n);
end

for j = 1:n-2
    for i = n:-1:j+2
        
        if anm
            subplot(1, 2, 1);
            heatmap(abs(A)); title('A'); getframe(gcf);
            subplot(1, 2, 2);
            heatmap(abs(B)); title('B'); getframe(gcf);
        end    
        
        [c, s] = givens(A(i-1,j), A(i,j));
        G = [c, s; -conj(s), conj(c)];
        A(i-1:i, j:n) = G * A(i-1:i, j:n);
        B(i-1:i, i-1:n) = G * B(i-1:i, i-1:n);
        if nargout >= 3
            Q(1:n, i-1:i) = Q(1:n, i-1:i) * G';
        end
        
        if anm
            subplot(1, 2, 1);
            heatmap(abs(A)); title('A'); getframe(gcf);
            subplot(1, 2, 2);
            heatmap(abs(B)); title('B'); getframe(gcf);
        end

        [c, s] = givens(-B(i, i), B(i, i-1));
        G = [c, s; -conj(s), conj(c)];
        B(1:i, i-1:i) = B(1:i, i-1:i) * G';
        B(i, i-1) = 0;
        A(1:n, i-1:i) = A(1:n, i-1:i) * G';
        if nargout == 4
            Z(1:n, i-1:i) = Z(1:n, i-1:i) * G';
        end
    end
    A(j+2:n, j) = 0;
end

if anm
    subplot(1, 2, 1);
    heatmap(abs(A)); getframe(gcf);
    subplot(1, 2, 2);
    heatmap(abs(B)); getframe(gcf);
end


