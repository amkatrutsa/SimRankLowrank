# SimRankLowrank

## About

This repository presents the approach to compute SimRank lowrank approximation [[I. Oseledets, G. Ovchinnikov, A. Katrutsa 2017]](http://comnet.oxfordjournals.org/content/early/2016/05/27/comnet.cnw008.abstract). SimRank is a similarity measure between graph vertices originally introduced in [[Jeh and Widom, 2002]](http://www-cs-students.stanford.edu/~glenj/simrank.pdf). But for large graphs storing SimRank for all pairs of vertices is intractable. Therefore we propose a method to store not exhaustive SimRank matrix, but only two lowrank matrices. Moreover, we provide an equation how to recover SimRank approximation from the adjacency matrix and two lowrank matrices.   

## Using in MATLAB

To test SimRankLowrank approximation approach, run the following commands in the directory you choose:
```
git clone https://github.com/amkatrutsa/SimRankLowrank.git
cd ./mcode
```
Then in MATLAB command line from the directory `mcode` run the command:
```
main
```

After that, you will have precise SimRank matrix `S`, computed by the naive SimRank algorithm, and its lowrank approximation `S_lr`, computed by our method.

## Using in C++

To use our code from C++ you need to install [Armadillo](http://arma.sourceforge.net/) library.
To test SimRankLowrank approximation approach,  run the following commands in the directory you choose:
```
git clone https://github.com/amkatrutsa/SimRankLowrank.git
cd ./ccode
mkdir build
cd ./build
cmake ..
make
```
After that you will have the binary file ``simrank_lowrank``.

## Test graphs

We carried out our experiments on the graphs from [DIMACS10 Collection](http://www.cise.ufl.edu/research/sparse/matrices/DIMACS10/index.html), which are in the directory `data`.
