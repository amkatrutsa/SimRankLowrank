PROJECT("simrank_lowrank")
CMAKE_MINIMUM_REQUIRED(VERSION 2.6.2)
SET(CMAKE_CXX_COMPILER /usr/bin/clang++)
LIST(APPEND CMAKE_CXX_FLAGS " -std=c++11 -Wall -O3") 
INCLUDE_DIRECTORIES(./src/)
ADD_EXECUTABLE(simrank_lowrank
main.cpp
./src/Parameters.cpp
./src/SimRankLowrank.cpp
)

TARGET_LINK_LIBRARIES(simrank_lowrank -larmadillo -lhdf5)
