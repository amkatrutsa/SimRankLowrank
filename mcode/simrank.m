function [S] = simrank(W, c, iter)
% This function implements the naive SimRank computation
% 
% Input:
% W - [n, n] - column normalized adjacency matrix
% c - double - constant from canonical SimRank
% iter - int - number of iteration
%
% Output:
% S - [n, n] - precise SimRank matrix
%
% Author: Alexandr Katrutsa, Skolkovo Institute of Science and Technology
% E-mail: aleksandr.katrutsa@phystech.edu

S = eye(size(W));
for k = 1:iter
    S = c * W' * S * W;
    S = S - diag(diag(S)) + eye(size(W));
end
end
