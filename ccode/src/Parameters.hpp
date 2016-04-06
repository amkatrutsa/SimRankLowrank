#include <string>

struct Parameters {
    Parameters();
    bool init(int argc, char* argv[]);
    void exit_with_help();
    // Filename with graph writing in MatrixMarket format, non-symmetrized
    std::string graph_filename;
    // SimRank parameter c
    double c;
    // Rank of approximation
    int rank;
    // Number of iteration in lowrank approximation loop
    int num_iter;
    // Oversampling parameter for probabilistic spectral decomposition
    int oversampling;
    // Filename for matrix U
    std::string U_filename;
    // Filename for matrix D
    std::string D_filename;
};
