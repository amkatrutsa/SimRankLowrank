function [S] = simrank_lowrank(W, c, r, p, k)
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
A1 = Wsq' * Wsq;
A2 = W' * diag(diag(W' * W)) * W;
for i = 1:k
    A3 = W' * U * D * U' * W;
    A = A1 - A2 + A3;
    M = A - diag(diag(A));
    [U, D] = ProbabSpectralDecomp(M, r, p, n);
end
S = eye(n) + W' * W - diag(diag(W' * W)) + U * D * U';
end

function [U, D] = ProbabSpectralDecomp(M, r, p, n)
Z = randn(n, r + p);
Y = M * Z;
Q = orth(Y);
B = Q' * (M * Q);
[V, Lambda, ~] = svd(B, 'econ');
U = Q * V;
D = Lambda;
end
