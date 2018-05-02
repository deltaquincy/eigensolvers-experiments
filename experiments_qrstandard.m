%% Experiments on Single-shift QR Algorithm for Schur Decomposition
%
% This file includes some experiments on Single-shift QR algorithm for Schur
% decomposition. Some of them performs comparisons on different selection on
% algorithm parameters, and the others shows the properties of the algorithm
% itself.
%
%  -------------------------------------------------
%  Experiments on Matrix Computations -- Spring 2018
%  Author: Liang Zilong
%  Date:   2018-04-01
%  -------------------------------------------------


%% Before experiments: Introduction
%
% This part introduces some special matrices in the following experiments.

clear; close all; clc;


% Block-tridiagonal matrix, introduced from five-point finite difference method
% for Poisson equations. It is symmetric and diagonally-dominant.
B = blocktridiag(3);

% Wilkinson's matrix. It is a built-in matrix in MATLAB, which is symmetric and
% tridiagonal. Some of its eigenvalues could be very close, but not equal.
% (see: https://blogs.mathworks.com/cleve/2013/04/15/wilkinsons-matrices-2/)
open wilkinson;
W21 = wilkinson(21);
plot(eig(W21), '.', 'MarkerSize', 20); grid on;
xlabel('i-th eigenvalue');
ylabel('\lambda');
title('Eigenvalues of Wilkinson matrix');


%% Experiment 1: Upper Hessenberg reduction using Householder or Givens
%
% This experiment performs a comparison between Householder transformation
% and Givens rotation on upper Hessenberg decomposition process.

clear; close all; clc;


n_max = 100;
repeat = 10;
matrices = {'general', 'block-tridiagonal', 'sparse', 'tridiagonal'};


for k = 1:4
    times_house = zeros(n_max, 1);
    times_givens = zeros(n_max, 1);
    
    for n = 1:n_max
        for j = 1:repeat
            switch k
                case 1  % General matrix (dense)
                    A = randn(n, n);
                case 2  % Block-tridiagonal
                    A = blocktridiag(n / 3);
                case 3  % Random, but sparse
                    A = triu(randn(n, n));
                    A = A + diag(randn(n-1, 1).*binornd(1, 0.4, [n-1, 1]), -1);
                    A = A + diag(randn(n-2, 1).*binornd(1, 0.4, [n-2, 1]), -2);
                case 4  % Wilkinson's
                    A = wilkinson(n);
            end
            
            t0_householder = clock;
            qrstandard_reduction(A);
            times_house(n) = times_house(n) + etime(clock, t0_householder);

            t0_givens = clock;
            hessenberg_givens(A);
            times_givens(n) = times_givens(n) + etime(clock, t0_givens);
        end
        
        times_house(n) = times_house(n) / repeat;
        times_givens(n) = times_givens(n) / repeat;
    end
    
    xx = zeros(n_max-9, 2);
    xx(:, 1) = log(10:n_max);
    xx(:, 2) = 1;
    yy_house = log(times_house(10:end));
    yy_givens = log(times_givens(10:end));
    coef_house = xx \ yy_house;
    coef_givens = xx \ yy_givens;
    
    figure(k);
    loglog(times_house, '.-'); grid on; hold on;
    loglog(times_givens, '.-');
    loglog(exp(xx(:, 1)), exp(coef_house(2) + coef_house(1) .* xx(:, 1)));
    text(0.7*exp(xx(1, 1)), exp(coef_house(2) + coef_house(1) .* xx(1, 1)), ...
         sprintf('k_{1} = %.3f', coef_house(1)));
    loglog(exp(xx(:, 1)), exp(coef_givens(2) + coef_givens(1) .* xx(:, 1)));
    text(0.7*exp(xx(1, 1)), exp(coef_givens(2) + coef_givens(1) .* xx(1, 1)), ...
         sprintf('k_{2} = %.3f', coef_givens(1)));
    xlabel('n');
    ylabel('time (s)');
    legend('Householder', 'Givens');
    title(sprintf('Experiment : Upper Hessenberg reduction (%s)', ...
          matrices{k}));
end


%% Experiment 2: Eigenvalues checking
%
% This experiment simply checks the wrote-up function 'qrstandard.m' works well.

clear; close all; clc;


n = 25;
matrices = {'real', 'real-symmetric', 'complex'};

for k = 1:3
    switch k
        case 1  % Random real matrix (normal)
            A = randn(n);
            xlabelstr = 'Re(\lambda)';
            ylabelstr = 'Im(\lambda)';
        case 2  % Wilkinson's matrix
            A = wilkinson(n);
            xlabelstr = 'i-th eigenvalue';
            ylabelstr = '\lambda';
        case 3  % Random complex matrix (normal)
            A = randn(n) + 1i*randn(n);
            xlabelstr = 'Re(\lambda)';
            ylabelstr = 'Im(\lambda)';
    end
    
    T = qrstandard(A);
    eigs_qrschur = diag(T);
    eigs_correct = eig(A);
    err = norm(sort(abs(eigs_qrschur)) - sort(abs(eigs_correct)), Inf);
    
    % Judge real or complex, and sort if real
    if isreal(eigs_correct)
        eigs_qrschur = sort(real(eigs_qrschur));
        eigs_correct = sort(eigs_correct);
        xlabelstr = 'i-th eigenvalue';
        ylabelstr = '\lambda';
    end
    
    figure(k);
    plot(eigs_qrschur, '.'); hold on;
    plot(eigs_correct, 'o');
    xlabel(xlabelstr);
    ylabel(ylabelstr);
    legend(sprintf('%s (qrschur)', matrices{k}), ...
           sprintf('%s (correct)', matrices{k}));
    title(sprintf('Experiment : Eigenvalues check (err: %.3e)', err));
end


%% Experiment 3: Comparison among different types of shift
%
% This experiment performs a comparison among different shifts on single-shifted
% QR algorithm. We choose non-shifted method, Rayleigh quotient shift and
% Wilkinson's shift to perform the comparison, and A is a random complex matrix.
%
% Caution: Performance of this experiment could be extremely low, thanks to the
%          non-shifted method.

clear; close all; clc;


n_max = 50;
repeat = 3;
shifts = {'none', 'rayleigh', 'wilkinson'};
legendcmd = 'legend(';

for k = 1:3
    times = zeros(n_max, 1);
    for n = 1:n_max
        for j = 1:repeat
            sprintf('k = %d, n = %d, repeat = %d', k, n, j)
            A = randn(n) + 1i * randn(n);
            
            t0 = clock;
            qrstandard(A, shifts{k}, 1e-5);
            times(n) = times(n) + etime(clock, t0);
        end
        times(n) = times(n) / repeat;
    end
    
    xx = zeros(n_max-1, 2);
    xx(:, 1) = log(2:n_max);
    xx(:, 2) = 1;
    yy = log(times(2:end));
    coef = xx \ yy;
    
    loglog(2:n_max, times(2:end), ...
           '.-', 'MarkerSize', 13, 'LineWidth', 1.5); grid on; hold on;
    loglog(exp(xx(:, 1)), exp(coef(2) + coef(1) .* xx(:, 1)));
    text(0.7*exp(xx(1, 1)), exp(coef(2) + coef(1) .* xx(1, 1)), ...
         sprintf('k_{%d} = %.3f', k, coef(1)));
    legendcmd = strcat(legendcmd, sprintf('''%s'',''%s (fit)'',', ...
                                          shifts{k}, shifts{k}));
end
xlabel('dim');
ylabel('average time');
legendcmd = strcat(legendcmd(1:end-1), ');');
eval(legendcmd);
title('Experiment : Comparison among different types of shift');


%% Experiment 4: Average iteration steps
%
% This experiment performs an observation on average iterations steps of single-
% shifted QR algorithm. Also, we choose three kind of matrices (random, block-
% tridiagonal and Wilkinson) to compare, to show the relationships between
% iterations steps and distributions of eigenvalues.

clear; close all; clc;


n = 25;
xx = 1:n;

matrices = {'rand', 'randn', 'blocktridiag', 'wilkinson'};

for k = 1:4
    switch k
        case 1  % Random matrix (uniform)
            A = rand(n, n);
        case 2  % Random matrix (normal)
            A = randn(n, n);
        case 3  % Block-tridiagonal matrix
            A = blocktridiag(n / 3);
        case 4  % Wilkinson's matrix
            A = wilkinson(n);
    end
    
    [eigs_qrschur, iters] = qrschur_observation(A);
    steps = iters - [iters(2:end); 0];
    avg = mean(steps);
    
    figure(k);
    [AX, H1, H2] = plotyy(xx, steps(end:-1:1), ...
                          xx, abs(eigs_qrschur(end:-1:1))); grid on; hold on;
    plot([0, n], [avg, avg], 'b-', 'LineWidth', 1.5);
    text(0.5, 1.05*avg, sprintf('AVG: %.3f', avg), 'FontWeight', 'bold');
    set(H1, 'Marker', '.', 'MarkerSize', 10, 'LineStyle', '--', 'Color', 'b');
    set(H2, 'Marker', '.', 'MarkerSize', 15, 'LineStyle', 'none');
    set(get(AX(1), 'ylabel'), 'string', 'iteration steps');
    set(get(AX(2), 'ylabel'), 'string', 'abs(\lambda)');
    xlabel('i-th eigenvalue');
    title(sprintf('Experiment: Average iteration steps (%s)', ...
          matrices{k}));
end


%% Appendix: Helper functions
%
% This part includes some helper functions for these experiments. Some of them
% are adapted from practical algorithms, some features added to show special
% properties of the algorithms.

function [H, U] = hessenberg_givens(A)
% -------------------------------------------------
% Helper function: upper Hessenberg decomposition 
%                  using Givens rotations
% -------------------------------------------------
    n = size(A);
    if n(1) ~= n(2)
        error('Fatal Error: The input is not a sqaure matrix!');
    end
    n = n(1);

    H = A;
    if nargout == 2
        U = eye(n, n);
    end

    for j = 1:n-2
        for i = n:-1:j+2
            if H(i, j) == 0 || abs(H(i, j)) < 2 * eps
                continue
            end
            [c, s, eta] = givens(H(i-1, j), H(i, j));
            G = [c, s; -conj(s), conj(c)];
            H(i-1:i, j) = [eta; 0];
            H(i-1:i, j+1:n) = G * H(i-1:i, j+1:n);
            H(:, i-1:i) = H(:, i-1:i) * G';
            if nargout == 2
                U(:, i-1:i) = U(:, i-1:i) * G';
            end
        end
    end
end

function [eigs_qrschur, iters] = qrschur_observation(A)
% -------------------------------------------------
% Helper function: Single-shifted QR algorithm
%                  with observation of convergense
%                  properties.
%                  (Adapted from 'qrstandard.m')
% -------------------------------------------------
    n = size(A);
    if n(1) ~= n(2)
        error('Fatal Error: The input is not a sqaure matrix!');
    end
    n = n(1);
    
    eigs_qrschur = zeros(n, 1);
    iters = zeros(n, 1);

    % Step 1: Upper Hessenberg reduction
    H = qrstandard_reduction(A);
    
    l = 0;
    m = 0;
    m_next = 0;
    iter = 0;
    while true
        % Step 2: Judge convergence of an eigenvalue
        for k = 1:n-1
            if abs(H(k+1, k)) <= eps * (abs(H(k, k)) + abs(H(k+1, k+1)))
                H(k+1, k) = 0;
            end
        end
        for m_next = m:n-2
            if H(n-m_next, n-m_next-1) ~=0
                break
            end
        end
        if m_next == n-2 && H(n-m, n-m-1) == 0
            m_next = n;
        end
        eigs_qrschur(n-m_next+1:n-m) = diag(H(n-m_next+1:n-m, n-m_next+1:n-m));
        iters(n-m_next+1:n-m) = iter;
        if m_next == n  % Final condition
            break
        end
        m = m_next;

        % Step 3: QR iteration
        [H22, G] = qriteration(H(l+1:n-m, l+1:n-m));
        H(l+1:n-m, l+1:n-m) = H22;
        for k = l+1:n-m-1
            H(1:l, k:k+1) = H(1:l, k:k+1) * G(:, :, k-l)';
            H(k:k+1, n-m+1:n) = G(:, :, k-l) * H(k:k+1, n-m+1:n);
        end
        iter = iter + 1;
    end
end

function A = blocktridiag(n)
% -------------------------------------------------
% Helper function: Generate 3*n-D block-tridiagonal
%                  matrix of five-point finite
%                  difference method.
% -------------------------------------------------
    A = - diag(ones(3*n-3,1), 3) - diag(ones(3*n-3,1), -3);
    D = [4, -1, 0; -1, 4, -1; 0, -1, 4];
    for k = 0:n-1
        A(3*k+1:3*k+3, 3*k+1:3*k+3) = D;
    end
end