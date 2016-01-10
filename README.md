# SimRankLowrank

## About

This repository presents the approach to compute SimRank lowrank approximation. SimRank is a measure of similarity between graph vertices originally introduced in [Jeh and Widom, 2002] (http://www-cs-students.stanford.edu/~glenj/simrank.pdf). 

## Running

To test SimRank Lowrank approximation approach, run the following commands in directory you choose:
```
git clone https://github.com/amkatrutsa/SimRankLowrank.git
cd ./mcode
```
Then in MATLAB command line from the directory `mcode` run the command:
```
main
```

After that you will have precise SimRank matrix `S` and its lowrank approximation `S_lr`.
