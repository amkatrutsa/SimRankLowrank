clear;
path_to_data = '../data/';
filename = strcat(path_to_data, 'delaunay_n10.mat');
W = load_data(filename);
[n, ~] = size(W);
% SimRank scale parameter, 0 < c < 1
c = 0.3;
% Overpampling parameter for Probabilistic Spectral Decomposition
p = 10;
% Number of iteration
k = 100;
% Rank approximation
r = 100;
W = W' + W; 
W = norm_by_col(W);
% Naive SimRank computation
S = simrank(W, c, k);
% Lowrank SimRank approximation 
S_lr = simrank_lowrank(W, c, r, p, k);