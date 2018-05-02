%% Experiments on QZ Algorithm for Periodic Schur Decomposition
%
% This file includes some experiments on QZ algorithm for periodic Schur
% Decomposition.
%
%  -------------------------------------------------
%  Experiments on Matrix Computations -- Spring 2018
%  Author: Liang Zilong
%  Date:   2018-04-23
%  -------------------------------------------------

%% Experiment 01: Time usage of periodic QZ algorithm
%
% This experiment checks out the time complexity of periodic QZ algorithm.
% Theoretically, it is O(p*n^3). So here we pick some 'p' values and 'n' values
% to test the time usage.

clear; close all; clc;

% Hyper-parameters
prange = [3, 10, 30];
plen = length(prange);
nmin = 3;
nmax = 100;
rep = 5;

% Legend parameters
lgdstr = cell(1, plen*2);

% Calculating
for i = 1:plen
    p = prange(i);
    t = zeros(nmax-nmin+1, 1);
    for n = nmin:nmax
        for k = 1:rep
            E = cell(1, p);
            A = cell(1, p);
            for j = 1:p
                E{j} = rand(n);
                A{j} = rand(n);
            end
            c = clock;
            qzperiodic(E, A);
            t(n-nmin+1) = t(n-nmin+1) + etime(clock, c);
        end
        t(n-nmin+1) = t(n-nmin+1) / rep;
    end
    
    % Semilog graph of time used
    loglog(nmin:nmax, t, '.-', 'LineWidth', 1.5); grid on; hold on;
    
    % Linear fit
    x = zeros(nmax-nmin+1, 2);
    x(:, 1) = log(nmin:nmax);
    x(:, 2) = 1;
    y = log(t);
    coef = x \ y;
    xx = exp(x(:, 1));
    yy = exp(coef(2) + coef(1) .* x(:, 1));
    loglog(xx, yy, '--', 'LineWidth', 1);
    text(xx(1), yy(1), sprintf('$k_{%d} = %.3f$', p, coef(1)), ...
         'Interpreter', 'latex', 'FontSize', 14);
    
    % Legend strings
    lgdstr{i*2-1} = sprintf('$p = %d$', p);
    lgdstr{i*2} = sprintf('$p = %d$ (linear fit)', p);
end

% Plotting
xlabel('dim');
ylabel('Time used (s)');
title('Experiment: Time usage of periodic $QZ$ algorithm');
legend(lgdstr);
setstyle(gca, 'latex');


%% Experiment 02: Time usage among different period

% Hyper-parameters
prange = 5:100;
plen = length(prange);
n = 30;
rep = 3;

% Placeholder
t = zeros(plen, 1);

% Calculating
for i = 1:plen
    p = prange(i);
    E = cell(1, p);
    A = cell(1, p);
    for k = 1:rep
        for j = 1:p
            E{j} = rand(n, n);
            A{j} = rand(n, n);
        end
        c = clock;
        qzperiodic(E, A);
        t(i) = t(i) + etime(clock, c);
    end
    t(i) = t(i) / rep;
end

% Linear fit
x = zeros(plen, 2);
x(:, 1) = prange;
x(:, 2) = 1;
y = t;
coef = x \ y;
xx = x(:, 1);
yy = coef(2) + coef(1) .* x(:, 1);

% Plotting
plot(prange, t, '.-', 'LineWidth', 1.2); grid on; hold on;
plot(xx, yy, '--', 'LineWidth', 0.8);
text(xx(1), yy(1)*1.3, sprintf('$k = %.3f$', coef(1)), ...
     'Interpreter', 'latex', 'FontSize', 14);
xlabel('Period');
ylabel('Time used');
title('Experiment: Time usage among different period');
setstyle(gca, 'latex');


%% Experiment 03: Combination shift strategy
%
% This experiment tests the combination shift strategy. For matrix sequences
% with different percentage of real eigenvalues, it compares the time when only
% using the double-shift and combination shift strategy.

clear; close all; clc;

% Hyper-parameters
p = 4;
n = 40;
rep = 4;
propslen = 30;
props = linspace(0, 0.98, propslen);

% Placeholders
E = cell(1, p);
A = cell(1, p);

% Time statistics
td = zeros(propslen, 1);
tc = zeros(propslen, 1);

for i = 1:propslen
    prop = props(i);
    for k = 1:rep
        for j = 1:p
            if j == 1
                nr = floor(n*prop);  % n-real
                nc = n - nr;         % n-complex
                A{1} = eigmat([rand(1, nr), rand(1, nc) + 1i*rand(1, nc)]);
            else
                A{j} = triu(rand(n, n));
            end
            E{j} = triu(rand(n, n));
        end
        c = clock;
        qzperiodic(E, A, 'double');
        td(i) = td(i) + etime(clock, c);

        c = clock;
        qzperiodic(E, A, 'comb');
        tc(i) = tc(i) + etime(clock, c);
    end
    td(i) = td(i) / rep;
    tc(i) = tc(i) / rep;
end

% Plotting
plot(props, tc./td, '+-', 'LineWidth', 1.5); grid on;
xlabel('Proportion of real eigenvalues');
ylabel('time(combination) / time(double)');
title('Experiment: Combination shift strategy');
setstyle(gca, 'latex');


%% Helper function

function A = eigmat(lambda)
% This is a helper function that generates a matrix with given eigenvalues.
    n = length(lambda);
    A = diag(lambda) + triu(rand(n, n), 1);
    Q = orth(randn(n, n));
    A = Q' * A * Q;
end