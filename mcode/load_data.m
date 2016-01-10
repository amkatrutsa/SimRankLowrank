function W = load_data(filename)
% Function loads data from mat-file downloaded from
% http://www.cise.ufl.edu/research/sparse/matrices/DIMACS10/index.html
% 
% Input
% filename - string - name of the file
% 
% Output
% W - [n,n] - adjacency matrix
% 
% Author: Alexandr Katrutsa, Skolkovo Institute of Science and Technology
% E-mail: aleksandr.katrutsa@phystech.edu
% Date: 20.11.2014

data = load(filename);
W = data.Problem.A;
end
