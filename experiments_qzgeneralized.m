%% Experiments on QZ Algorithm for Generalized Schur Decomposition
%
% This file includes some experiments on QZ algorithm for generalized Schur
% Decomposition.
%
%  -------------------------------------------------
%  Experiments on Matrix Computations -- Spring 2018
%  Author: Liang Zilong
%  Date:   2018-04-14
%  -------------------------------------------------


%% Animation demo: Hessenberg-triangular reduction
%
% This part shows the animation of the process of Hessenberg-triangular
% reduction.

clear; close all; clc;

A = rand(15);
B = rand(15);
qzgeneralized_reduction(A, B, true);


%% Animation demo: Single-shift QZ step
%
% This part shows the animation of the process of single-shift QZ step.

clear; close all; clc;

A = rand(15);
B = rand(15);
[A, B] = qzgeneralized_reduction(A, B);
qziteration_anm(A, B, 'single');


%% Animation demo: Double-shift QZ step
%
% This part shows the animation of the process of double-shift implicit QZ step.

clear; close all; clc;

A = rand(15);
B = rand(15);
[A, B] = qzgeneralized_reduction(A, B);
qziteration_anm(A, B, 'double');


%% Experiment: Single-shift v.s. double-shift implicit QZ
%
% Ward (1975) raised Combination Shift QZ algorithm. This experiment carried out
% a comparison between Single-shift, double-shift and combination-shift QZ
% algorithm.

clear; close all; clc;

n_max = 30;
times = zeros(3, n_max);
rep = 3;

for j = 1:3
    switch j
        case 1
            shift = 'single';
        case 2
            shift = 'double';
        case 3
            shift = 'comb';
    end
    for i = 1:n_max
        for k = 1:rep
            A = randn(i, i);
            B = eye(i);
            t = clock;
            qzgeneralized(A, B, shift);
            times(j, i) = times(j, i) + etime(clock, t);
        end
        times(j, i) = times(j, i) / rep;
    end
    semilogy(times(j, :), '.-', 'LineWidth', 1, 'MarkerSize', 10); hold on;
end
title('Experiment: Comparison between shift types');
xlabel('$\mathrm{dim}(A)$', 'Interpreter', 'latex');
ylabel('Time used to get quasi-triangular forms');
legend('single', 'double', 'combiniation');
setstyle(gca, 'latex');


%% Experiment: QZ algorithm without deflation
%
% Moler and Stewart (1973) raised an example that illustrates QZ algorithm is
% unaffected by an ill-conditioned B matrix. Here we adapt this example.

clear; close all; clc;

% The same matrix pair as (Moler and Stewart, 1973).
A = [50 -60 50 -27  6 6;
     38 -28 27 -17  5 5;
     27 -17 27 -17  5 5;
     27 -28 38 -17  5 5;
     27 -28 27 -17 16 5;
     27 -28 27 -17 5 16;];
B = [16 5 5 5 -6 5;
     5 16 5 5 -6 5;
     5 5 16 5 -6 5;
     5 5 5 16 -6 5;
     5 5 5 5 -6 16;
     6 6 6 6 -5  6;];

[A1, B1] = qzschur_withoutdeflation(A, B);
eig_qzschur = [eig(A1(1:2, 1:2), B1(1:2, 1:2));
               eig(A1(3:4, 3:4), B1(3:4, 3:4));
               eig(A1(5:6, 5:6), B1(5:6, 5:6));];
eig_qzschur = sort(eig_qzschur);
eig_right = sort(eig(A, B));
eig_err = norm(eig_right(1:4) - eig_qzschur(1:4), Inf);


%% Experiment: A further step, when B is singular
%
% This experiment carried out another example illustrating that QZ algorithm is
% unaffected by an singular B matrix.

clear; close all; clc;

A = rand(20);
B = rand(20); B(7, :) = 0;

[A1, B1] = qzgeneralized(A, B);
[A2, B2] = qzschur_withoutdeflation(A, B);
eig_qzschur = sort(eig(A1, B1));
eig_qzschur_withoutinf = eig_qzschur(1:end-1);
eig_qzschurwd = sort(eig(A2, B2));
eig_qzschurwd_withoutinf = eig_qzschurwd(1:end-1);
eig_err = norm(eig_qzschurwd_withoutinf - eig_qzschur_withoutinf, Inf);

%% Experiment: Polynomial eigenvalue problems
% 
% This experiment is an example on polynomial eigenvalue problems, which is an
% application of generalized eigenvalue problems.

clear; close all; clc;

n = 50;
a = @(x) x.^2 .* (pi - x).^2 - 2.7;

M = (pi/2) * eye(n);
K = (pi/2) * diag((1:n).^2);
C = zeros(n, n);
for i = 1:n
    for j = 1:n
        C(i, j) = 0.1 * integral(@(x) a(x).*sin(i.*x).*sin(j.*x), 0, pi);
    end
end
[A, B] = pepmat(K, C, M);
[At, Bt] = qzgeneralized(A, B);
eigs1 = eig(At, Bt);
eigs2 = polyeig(K, C, M);
err = norm(sort(abs(eigs1))-sort(abs(eigs2)));

plot(eigs1, '.'); hold on; 
plot(eigs2, 'o');
xlabel('Re($\lambda$)', 'Interpreter', 'latex');
ylabel('Im($\lambda$)', 'Interpreter', 'latex');
text(-0.05, 0, sprintf('err = $%.3e$', err), 'Interpreter', 'latex', 'FontSize', 15);
title(sprintf('Experiment: PEP on wave equation ($n = %d$)', n));
legend('QZ', 'polyeig');
setstyle(gca, 'latex');


%% Appendix: Helper functions

function A = eigenfixed(lambda)
    n = length(lambda);
    A = diag(lambda) + triu(randn(n, n), 2);
    Q = orth(randn(n, n));
    A = Q' * A * Q;
end

function [A, B] = pepmat(varargin)
    m = length(varargin) - 1;
    n = length(varargin{1});
    A = diag(-ones((m-1)*n, 1), -n);
    B = eye(m*n);
    for i = 1:m
        A((i-1)*n+1:i*n, (m-1)*n+1:m*n) = varargin{i};
    end
    B((m-1)*n+1:m*n, (m-1)*n+1:m*n) = varargin{m+1};
    B = -B;
end

function [A, B, Q, Z] = qzschur_withoutdeflation(A, B, tol)
% -------------------------------------------------
% Helper function: QZ algorithm without the process
%                  of deflation.
%                  (Adapted from qzschur.m)
% -------------------------------------------------
    if nargin < 3
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
        [A, B] = htreduction(A, B);
    elseif nargout == 3
        [A, B, Q] = htreduction(A, B);
    elseif nargout == 4
        [A, B, Q, Z] = htreduction(A, B);
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
        if m == n-3 && (A(2,1) == 0 || A(3, 2) == 0)  % Final condition
            break
        end
        l = 0;
        for k = n-m-2:-1:1
            if A(k+1,k) == 0
                l = k;
                break
            end
        end

        % Step 3: QZ iteration
        if nargout <= 2
            [A, B] = qziteration(A, B, l, m, 'comb');
        elseif nargout == 3
            [A, B, Q] = qziteration(A, B, l, m, Q, 'comb');
        elseif nargout == 4
            [A, B, Q, Z] = qziteration(A, B, l, m, Q, Z, 'comb');
        end
    end
end

function [A, B] = qziteration_anm(A, B, shift)
% -------------------------------------------------
% Helper function: QZ iteration step with animation
%                  (Adapted from qziteration.m)
% -------------------------------------------------

    n = length(A);
    l = 0;
    m = 0;
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
            subplot(1, 2, 1);
            heatmap(abs(A)); title('A'); getframe(gcf);
            subplot(1, 2, 2);
            heatmap(abs(B)); title('B'); getframe(gcf);

            [c, s] = givens(A(k, k), A(k+1, k));
            G = [c, s; -conj(s), conj(c)];
            A(k:k+1, k:n) = G * A(k:k+1, k:n);
            B(k:k+1, k:n) = G * B(k:k+1, k:n);
            if nargout >= 3
                Q(1:n, k:k+1) = Q(1:n, k:k+1) * G';
            end
            
            subplot(1, 2, 1);
            heatmap(abs(A)); title('A'); getframe(gcf);
            subplot(1, 2, 2);
            heatmap(abs(B)); title('B'); getframe(gcf);
            [c, s] = givens(-B(k+1, k+1), B(k+1, k));
            GZ(:, :, k-l) = [c, s; -conj(s), conj(c)];
            B(1:k+1, k:k+1) = B(1:k+1, k:k+1) * GZ(:, :, k-l)';
            B(k+1, k) = 0;
            if nargout == 4
                Z(1:n, k:k+1) = Z(1:n, k:k+1) * GZ(:, :, k-l)';
            end
        end
        for k = l+1:n2-1
            subplot(1, 2, 1);
            heatmap(abs(A)); title('A'); getframe(gcf);
            subplot(1, 2, 2);
            heatmap(abs(B)); title('B'); getframe(gcf);
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
            subplot(1, 2, 1);
            heatmap(abs(A)); title('A'); getframe(gcf);
            subplot(1, 2, 2);
            heatmap(abs(B)); title('B'); getframe(gcf);

            % Z-step
            w = house(B(k+2, k:k+2));  % TODO: è¿???????ä¸??¶ï?ä¸???è¦????¢æ??ä»¥ç?´æ?¥å½¢??äº? Householder ??
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
            subplot(1, 2, 1);
            heatmap(abs(A)); title('A'); getframe(gcf);
            subplot(1, 2, 2);
            heatmap(abs(B)); title('B'); getframe(gcf);

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
        subplot(1, 2, 1);
        heatmap(abs(A)); title('A'); getframe(gcf);
        subplot(1, 2, 2);
        heatmap(abs(B)); title('B'); getframe(gcf);

        [c, s] = givens(-B(n2, n2), B(n2, n2-1));
        G = [c, s; -conj(s), conj(c)];
        A(1:n2, n2-1:n2) = A(1:n2, n2-1:n2) * G';
        B(1:n2, n2-1:n2) = B(1:n2, n2-1:n2) * G';
        if nargout == 4
            Z(1:n, n2-1:n2) = Z(1:n, n2-1:n2) * G';
        end
        subplot(1, 2, 1);
        heatmap(abs(A)); title('A'); getframe(gcf);
        subplot(1, 2, 2);
        heatmap(abs(B)); title('B'); getframe(gcf);

        A = triu(A, -1);
        B = triu(B);

    end
end