function [w, gamma] = house(a)
% HOUSE    Compute Householder transformation
%
% A Householder matrix: H = I-ww' denotes a unitray transformation. This
% function computes the vector 'w' and a scalar 'gamma', given a vector 'a',
% so that Ha = gamma*e_1. Note that norm(w, 2) = sqrt(2).
%
% argin:
%   a - A column vector to transform (real or complex)
%
% argout:
%   w     - The vector in the Householder transformation H = I-ww'
%   gamma - The scalar in the formula Ha = gamma*e_1
%
% -------------------------------------------------
% Experiments on Matrix Computations -- Spring 2018
% Author: Liang Zilong
% Date:   2018-03-31
% -------------------------------------------------

if isrow(a)
    a = a';
end

w = a;
gamma = norm(a, 2);

if gamma == 0
    w(1) = sqrt(2);
    return
end

if w(1) ~= 0
    delta = conj(w(1)) / abs(w(1));
else
    delta = 1;
end

w = (delta / gamma) * w;
w(1) = w(1) + 1;

w = w / sqrt(w(1));
gamma = - conj(delta) * gamma;

