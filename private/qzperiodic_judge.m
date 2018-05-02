function [l, m] = qzperiodic_judge(E, A)
% QZPERIODIC_JUDGE    Judge convergence and perform block-deflation for periodic
%                     QZ algorithm.
% 
% This is a helper function that determines block-deflation parameters l and m
% for periodic QZ algorithm.
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-04-22
% -------------------------------------------------

% Check inputs
p = length(E);
n = length(E{1});

% Initialize
l = 0;
m = 0;

for k = n:-1:3
    if A{1}(k, k-1) ~= 0
        if A{1}(k-1, k-2) ~= 0
            m = n - k;
            break;
        else
            H = eye(2, 2);
            for j = p:-1:1
                H = H / E{j}(k-1:k, k-1:k);
                H = H * A{j}(k-1:k, k-1:k);
            end
            if isreal(eig(H))
                m = n - k;
                break;
            end
        end
    end
    m = n - k;
end

if m == n-3
    if A{1}(2, 1) == 0 && A{1}(3, 2) == 0
        m = n;
    elseif A{1}(2, 1) == 0
        H = eye(2, 2);
        for j = p:-1:1
            H = H / E{j}(2:3, 2:3);
            H = H * A{j}(2:3, 2:3);
        end
        if ~isreal(eig(H))
            m = n;
        end
    elseif A{1}(3, 2) == 0
%         m = n - 2;
%         H = eye(2, 2);
%         for j = p:-1:1
%             H = H / E{j}(1:2, 1:2);
%             H = H * A{j}(1:2, 1:2);
%         end
%         if ~isreal(eig(H))
%             m = n;
%         end
        m = n;
    end
end

for k = 1:n-m-1
    if A{1}(k+1, k) == 0
        l = k;
    end
end

