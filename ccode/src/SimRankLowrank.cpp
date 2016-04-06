#include "SimRankLowrank.hpp"
#include <armadillo>
#include <iostream>

SimRankLowrank::SimRankLowrank() {}

SimRankLowrank::SimRankLowrank(SpMat& scaled_adjacency_mat, int num_iter,
                                         double c, int rank, int oversampling_p) :
                                SimRankBase(c, num_iter), 
                                scaled_adjacency_mat_(scaled_adjacency_mat),
                                rank_(rank),
                                oversampling_(oversampling_p) {}

bool SimRankLowrank::init(const Parameters& p) {
    if ((p.c <= 0) && (p.c >= 1)) {
        printf("SimRank parameter c has to lie in segment (0, 1)\n");
        return false;
    }
    c_ = p.c;
    
    if (p.num_iter <= 0) {
        printf("Number of iterations has to be positive\n");
        return false;
    }
    num_iter_ = p.num_iter;
    graph_filename_ = p.graph_filename;
    U_filename_ = p.U_filename;
    D_filename_ = p.D_filename;
    
    if (!ReadGraph()){
        printf("Can not read graph form %s\n", graph_filename_.c_str());
        return false;
    }
    ScaleAdjacencyMatrix();
    
    if (p.oversampling < 0) {
        printf("Oversampling parameter has to be positive\n");
        return false;
    }
    if (p.oversampling >= num_vertex_) {
        printf("Oversampling parameter has to be less than the dimension of "
                "the adjacency matrix = %zu\n", num_vertex_);
        return false;
    }
    oversampling_ = p.oversampling;
    
    if (p.rank <= 0) {
        printf("Rank of approximation has to be positive\n");
        return false;
    }
    if (p.rank >= num_vertex_) {
        printf("Rank of approximation has to be less than the dimension of "
                "the adjacency matrix = %zu\n", num_vertex_);
        return false;
    }
    rank_ = p.rank;
    return true;
}

bool SimRankLowrank::ReadGraph() {
    printf("Reading graph from %s...\n", graph_filename_.c_str());
    std::ifstream graph_file(graph_filename_);
    std::string current_line; 
    if (!graph_file.is_open()) {
        std::cout << "Can not open file " << graph_filename_ << std::endl; 
        return false;
    }
    arma::umat loc;
    int i = 0;
    while(std::getline(graph_file, current_line)) {
        if (current_line[0] == '%')
            continue;
        std::istringstream iss(current_line);
        if (num_vertex_ == 0 && num_edges_ == 0) {
            iss >> num_vertex_ >> num_vertex_ >> num_edges_;
            loc.reshape(2, num_edges_);
        }
        else {
            int first_vertex, second_vertex;
            iss >> first_vertex >> second_vertex;
            loc(0, i) = first_vertex - 1;
            loc(1, i) = second_vertex - 1;
            ++i;
        }
    }
    Vec identity(num_edges_, arma::fill::ones);
    scaled_adjacency_mat_ = SpMat(loc, identity, num_vertex_, num_vertex_);
    scaled_adjacency_mat_ += scaled_adjacency_mat_.t();
    graph_file.close();
    printf("Number of vertices = %zu\n", num_vertex_);
    printf("Number of edges = %zu\n", num_edges_);
    printf("Reading graph from %s...Done\n", graph_filename_.c_str());
    return true;
}

void SimRankLowrank::ScaleAdjacencyMatrix() {
    Vec col_sum(num_vertex_, arma::fill::zeros);
    for (SpMat::const_iterator it = scaled_adjacency_mat_.begin(); it != scaled_adjacency_mat_.end(); ++it)
        col_sum(it.col()) += *it;
    for (SpMat::iterator it = scaled_adjacency_mat_.begin(); it != scaled_adjacency_mat_.end(); ++it)
        (*it) /= col_sum(it.col());
    scaled_adjacency_mat_ *= sqrt(c_);
}

bool SimRankLowrank::compute() {
    printf("Compute low rank approximation of SimRank...\n");
    U_ = Mat(num_vertex_, rank_ + oversampling_, arma::fill::randu);
    d_ = Vec(rank_ + oversampling_, arma::fill::ones);
    // Compute diag(A1) through many matrix by vector multiplications
    SpMat C_sq = scaled_adjacency_mat_ * scaled_adjacency_mat_;
    SpMat A1 = C_sq.t() * C_sq;
    // Compute diagonal of CTC
    SpMat CTC = scaled_adjacency_mat_.t() * scaled_adjacency_mat_;
    // Compute diag(A2) through many matrix by vector multiplications
    Vec diagonal_CTC = CTC.diag();
    arma::umat loc(2, scaled_adjacency_mat_.n_rows);
    for (size_t i = 0; i < loc.n_cols; ++i) {
        loc(0, i) = i;
        loc(1, i) = i;
    }
    SpMat diag_CTC(loc, diagonal_CTC, num_vertex_, num_vertex_);
    SpMat A2 = scaled_adjacency_mat_.t() * diag_CTC * scaled_adjacency_mat_;
    for (size_t i = 0; i < num_iter_; ++i) {
        std::cout << "Iteration " << i + 1 << std::endl;
        Mat UtW = U_.t() * scaled_adjacency_mat_;
        Mat A3 = UtW.t() * arma::diagmat(d_) * UtW;
        // Vec diag_A4 = Vec(A1.diag()) - Vec(A2.diag()) + Vec(A3.diag());
        Mat A = Mat(A1) - Mat(A2) + A3; 
        A -= arma::diagmat(A.diag());
        printf("SVD...\n");
        if (!ProbabilisticSpectralDecompositionSmall(A)) {
            printf("Error in probabilistic spectral decomposition\n");
            return false;
        }
        printf("SVD...Done\n");
    }
    printf("Compute low rank approximation of SimRank... Done\n");
    return true;
}

bool SimRankLowrank::save() {
    printf("Saving lowrank matrices...\n");
    if (!U_.save(U_filename_, arma::csv_ascii)) {
        printf("Can not save U matrix in %s\n", U_filename_.c_str());
        return false;
    }
    if (!d_.save(D_filename_, arma::csv_ascii)) {
        printf("Can not save D matrix in %s\n", D_filename_.c_str());
        return false;
    }
    printf("Saving lowrank matrices...Done\n");
    return true;
}


bool SimRankLowrank::ProbabilisticSpectralDecompositionSmall(const Mat& A) {
    // Generate matrix with elements from standard normal distribution
    int n = A.n_rows;
    Mat Z(n, rank_ + oversampling_, arma::fill::randn);
    Mat Y = A * Z;
//    Mat Q;
//    printf("Orth...\n");
//    if (!arma::orth(Q, Y)) {
//        printf("Error in orthonormal basis generation\n");
//        return false;
//    }
//    printf("Orth...Done\n");
//    printf("Dim Q = (%d, %d)\n", Q.n_rows, Q.n_cols);
//    Mat I = arma::eye(Q.n_cols, Q.n_cols);
//    Mat D = I - Q.t() * Q;
//    printf("D norm = %lf\n", arma::norm(D));
    Mat Q, R;
    printf("QR...\n");
    if (!arma::qr_econ(Q, R, Y)) {
        printf("Error in orthonormal basis generation\n");
        return false;
    }
    printf("QR...Done\n");
//    printf("Dim Q1 = (%d, %d)\n", Q1.n_rows, Q1.n_cols); 
//    D = I - Q1.t() * Q1;
//    printf("D norm = %lf\n", arma::norm(D));
//    D = Q - Q1;
//    printf("D norm = %lf\n", arma::norm(D));
    Mat B = Q.t() * (A * Q);
    Mat U1;
    printf("Eigen decomposition...\n");
    if (!arma::eig_sym(d_, U1, B)) {
        printf("Error in eigen decomposition\n");
        return false;
    }
    printf("Eigen decomposition...Done\n");
    U_ = Q * U1;
    return true;
}

//VectorXd SimRankLowrank::matvec(const VectorXd& v, const MatrixXd& U, const VectorXd& diag_D,
//                                const SpMat& A1, const SpMat& A2) {
//    int vector_dim = scaled_adjacency_mat_.rows();
//    // Compute A1 * v
//    VectorXd res1 = A1 * v; //matvec_A1(v);
//    // Compute A2 * v
//    VectorXd res2 = A2 * v; //matvec_A2(v, diagonal_CTC);
//    // Compute A3 * v
//    VectorXd res3 = matvec_A3(v, U, diag_D);
//    // Compute diag(A3)
//    VectorXd diagonal_A3(vector_dim);
//    MatrixXd A3_sqrt = diag_D.cwiseSqrt().asDiagonal() * U.transpose() * scaled_adjacency_mat_;
//    for (int i = 0; i < vector_dim; ++i)
//        diagonal_A3(i) = A3_sqrt.transpose().row(i) * A3_sqrt.col(i);
//    VectorXd res4 = A1.diagonal() - A2.diagonal() + diagonal_A3;
//    res4.cwiseProduct(v); 
//    return res1 - res2 + res3 - res4;
//}
//
//VectorXd SimRankLowrank::matvec_A1(const VectorXd& v) {
//    VectorXd res(scaled_adjacency_mat_.rows());
//    res = scaled_adjacency_mat_ * (scaled_adjacency_mat_ * v);
//    res = scaled_adjacency_mat_.transpose() * (scaled_adjacency_mat_.transpose() * res);
//    return res;
//}
//
//VectorXd SimRankLowrank::matvec_A2(const VectorXd& v, const VectorXd& diagonal_CTC) {
//    VectorXd res = scaled_adjacency_mat_ * v;
//    res = res.cwiseProduct(diagonal_CTC);
//    res = scaled_adjacency_mat_.transpose() * res;
//    return res;
//}
//
//VectorXd SimRankLowrank::matvec_A3(const VectorXd& v, const MatrixXd& U, const VectorXd& diag_D) {
//    VectorXd res = U.transpose() * (scaled_adjacency_mat_ * v);
//    res = res.cwiseProduct(diag_D);
//    res = scaled_adjacency_mat_.transpose() * (U * res);
//    return res;
//}
