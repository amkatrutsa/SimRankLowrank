#include "SimRankBase.hpp"
#include "Parameters.hpp"
#include <armadillo>

typedef arma::sp_mat SpMat;
typedef arma::mat Mat;
typedef arma::vec Vec;

// This class represents the SimRankLowrank approximation algorithm
class SimRankLowrank : public SimRankBase {
    public:
        SimRankLowrank();
        SimRankLowrank(SpMat& scaled_adjacency_mat, int num_iter, double c,
                            int rank, int oversampling_p);
        bool init(const Parameters& p);
        bool compute();
        bool save();
    private:
        virtual bool ReadGraph();
        void ScaleAdjacencyMatrix();
        void ProbabilSpectralDecomposition(Mat& U, Vec& diag_D, const SpMat& scaled_adj_matrix, 
                                           int rank, int oversample_p, const SpMat& A1, const SpMat& A2);
        bool ProbabilisticSpectralDecompositionSmall(const Mat& A);
        Vec matvec(const Vec& v, const Mat& U, const Vec& diag_D,
                        const SpMat& A1, const SpMat& A2);
        // Compute A1 * v
        Vec matvec_A1(const Vec& v);
        // Compute A2 * v
        Vec matvec_A2(const Vec& v, const Vec& diagonal_CTC);
        // Compute A3 * v
        Vec matvec_A3(const Vec& v, const Mat& U, const Vec& diag_D);
        SpMat scaled_adjacency_mat_;
        size_t rank_;
        size_t oversampling_;
        size_t num_iter_;
        // Lowrank rectangular matrix
        Mat U_;
        // Diagonal for the diagonal matrix
        Vec d_;
        std::string graph_filename_;
        std::string U_filename_;
        std::string D_filename_;
        size_t num_vertex_ = 0;
        size_t num_edges_ = 0;
};