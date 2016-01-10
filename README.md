# SimRankLowrank

## About

This repository presents the approach to compute SimRank lowrank approximation. SimRank is a similarity measure between graph vertices originally introduced in [[Jeh and Widom, 2002]] (http://www-cs-students.stanford.edu/~glenj/simrank.pdf). But for large graphs storing SimRank for all pairs of vertices is intractable. Therefore we propose a method to store not exhaustive SimRank matrix, but only two lowrank matrices. Moreover, we provide an equation how to recover SimRank approximation from the adjacency matrix and two lowrank matrices.   

## Running

To test SimRankLowrank approximation approach, run the following commands in the directory you choose:
```
git clone https://github.com/amkatrutsa/SimRankLowrank.git
cd ./mcode
```
Then in MATLAB command line from the directory `mcode` run the command:
```
main
```

After that, you will have precise SimRank matrix `S`, computed by naive SimRank algorithm, and its lowrank approximation `S_lr`, computed by our method.

## Test graphs

We carried out our experiment on the graphs from [DIMACS10 Collection](http://www.cise.ufl.edu/research/sparse/matrices/DIMACS10/index.html), which are in the directory `data`.
