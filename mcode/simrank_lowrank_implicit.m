function [S, U, D] = simrank_lowrank_implicit(W, c, r, p, k, n_iter_diag, ...
                                                return_S)
% Function computes SimRank base on lowrank approximation
%
% S = simrank_lowrank(W, c, r, p, k) returns a SimRank matrix approximation
%
% Input:
% W - [n, n] - column normalized adjacency matrix
% c - double - constant from canonical SimRank
% r - int - appropriate approximation rank
% p - int - oversampling parameter
% k - int - number of iteration
% n_iter_diag - int - number of iterations to estimate diagonal in matrix
% A
% return_S - bool - if True return Simrank matrix, otherwise return only
% factors U and D
%
% Output
% S - [n, n] - SimRank lowrank approximation
% or
% U [n, r+p] and D [r+p, r+p] - the matrices, which can reconstruct
% lowrank form of the matrix S
%
% Author: Alexandr Katrutsa, Skolkovo Institute of Science and Technology
% E-mail: aleksandr.katrutsa@phystech.edu
% Date: 19.11.2014
% 

W = sqrt(c) * W;
[n, ~] = size(W);
U = ones(n, r + p);
D = eye(r + p);
Wsq = W * W;
matvec_A1 = @(x) Wsq' * (Wsq * x);
matvec_A2 = @(x) W' * (diag(W' * W) .* (W * x));
for i = 1:k
    matvec_A3 = @(x) W' * (U * (diag(D) .* (U' * (W * x))));

    matvec_A = @(x) matvec_A1(x) - matvec_A2(x) + matvec_A3(x);
    A_diag = stoch_diag_est(matvec_A, n, n_iter_diag);
    matvec = @(x) matvec_A1(x) - matvec_A2(x) + matvec_A3(x) - A_diag .* x;
    [U, D] = ProbabSpectralDecomp_memeff(matvec, r, p, n);
end
if return_S
    S = eye(n) + W' * W - diag(diag(W' * W)) + U * D * U';
else
    S = NaN;
end
end

function [d] = stoch_diag_est(matvec, n, num_iter)
t = zeros(n, 1);
q = zeros(n, 1);
for i = 1:num_iter
    v = 2 * randi([0, 1], n, 1) - 1;
    t = t + v .* matvec(v);
    q = q + v .* v; 
end
d = t ./ q;
end

function [U, D] = ProbabSpectralDecomp_memeff(matvec, r, p, n)
Z = randn(n, r + p);
Y = matvec(Z);
Q = orth(Y);
B = Q' * (matvec(Q));
[V, Lambda, ~] = svd(B, 'econ');
U = Q * V;
D = Lambda;
end

