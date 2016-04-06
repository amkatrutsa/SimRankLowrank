#include "src/SimRankLowrank.hpp"

int main(int argc, char* argv[]) {
    arma::wall_clock t;
    t.tic();
    Parameters p;
    if (!p.init(argc, argv))
        return EXIT_FAILURE;
    SimRankLowrank s;
    if (!s.init(p))
        return EXIT_FAILURE;
    if (!s.compute())
        return EXIT_FAILURE;
    if (!s.save())
        return EXIT_FAILURE;
    double time = t.toc();
    printf("Elapsed time = %f\n", time);
    return EXIT_SUCCESS;
}