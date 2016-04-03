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
    
    if (p.oversampling < 0) {
        printf("Oversampling parameter has to be positive\n");
        return false;
    }
    oversampling_ = p.oversampling;
    
    if (p.rank <= 0) {
        printf("Rank of approximation has to be positive\n");
        return false;
    }
    rank_ = p.rank;
    
    if (!ReadGraph()){
        printf("Can not read graph form %s\n", graph_filename_.c_str());
        return false;
    }
    ScaleAdjacencyMatrix();
    return true;
}

bool SimRankLowrank::ReadGraph() {
    return true;
}

void SimRankLowrank::ScaleAdjacencyMatrix() {
    
}

bool SimRankLowrank::compute() {
    compute_small();
    return true;
}

void SimRankLowrank::compute_small() {
    printf("Compute low rank approximation of SimRank...\n");
    int num_vertex = scaled_adjacency_mat_.n_cols;
    Mat U(num_vertex, rank_ + oversampling_, arma::fill::randu);
    Vec diagonal_vec(rank_ + oversampling_, arma::fill::ones);
    // Compute diag(A1) through many matrix by vector multiplications
    SpMat A1 = scaled_adjacency_mat_.t() * scaled_adjacency_mat_.t() * 
               scaled_adjacency_mat_ * scaled_adjacency_mat_;
    // Compute diagonal of CTC
    SpMat CTC = scaled_adjacency_mat_.t() * scaled_adjacency_mat_;
    // Compute diag(A2) through many matrix by vector multiplications
    Vec diagonal_CTC = CTC.diag();
    arma::umat loc(scaled_adjacency_mat_.n_rows, 2);
    for (size_t i = 0; i < loc.n_rows; ++i) {
        loc(i, 0) = i;
        loc(i, 1) = i;
    }
    SpMat diag_CTC(loc, diagonal_CTC, scaled_adjacency_mat_.n_rows, scaled_adjacency_mat_.n_rows);
    SpMat A2 = scaled_adjacency_mat_.t() * diag_CTC * scaled_adjacency_mat_;
    for (size_t i = 0; i < num_iter_; ++i) {
        std::cout << "Iteration " << i + 1 << std::endl;
        Mat A3 = scaled_adjacency_mat_.t() * U * arma::diagmat(diagonal_vec) * U.t() * scaled_adjacency_mat_;
        Vec diag_A4 = Vec(A1.diag()) - Vec(A2.diag()) + Vec(A3.diag());
        Mat A = Mat(A1) - Mat(A2) + A3; 
        A -= arma::diagmat(diag_A4);
        Mat V;
        arma::svd_econ(U, diagonal_vec, V, A);
    }
//    Vec identity_vec(num_vertex, arma::fill::ones);
//    SpMat scaled_adj_mat_product = scaled_adjacency_mat_.t() * scaled_adjacency_mat_;
//    simrank_ = Mat(scaled_adj_mat_product);
//    simrank_.diag() = identity_vec;
//    simrank_ += U * arma::diagmat(diagonal_vec) * U.t();
    printf("Compute low rank approximation of SimRank... Done\n");
}

bool SimRankLowrank::save() {
    return true;
}

//void SimRankLowrank::compute_large() {
//    printf("Compute low rank approximation of SimRank...\n");
//    int num_vertex = scaled_adjacency_mat_.cols();
//    MatrixXd U = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Random(num_vertex, rank_ + oversampling_p_);
//    VectorXd diagonal_vec(rank_ + oversampling_p_);
//    for (int i = 0; i < rank_ + oversampling_p_; ++i)
//        diagonal_vec(i) = 1;
//    // Compute diag(A1) through many matrix by vector multiplications
//    SpMat A1 = scaled_adjacency_mat_.transpose() * scaled_adjacency_mat_.transpose() * 
//               scaled_adjacency_mat_ * scaled_adjacency_mat_;
//    VectorXd diagonal_A1 = A1.diagonal();
//    // Compute diagonal of CTC
//    SpMat CTC = scaled_adjacency_mat_.transpose() * scaled_adjacency_mat_;
//    // Compute diag(A2) through many matrix by vector multiplications
//    VectorXd diagonal_CTC = CTC.diagonal();
//    SpMat diag_CTC(num_vertex, num_vertex);
//    std::vector<T> triplet_list;
//    triplet_list.reserve(num_vertex);
//    for (int i = 0; i < num_vertex; ++i)
//        triplet_list.push_back(T(i, i, diagonal_CTC(i)));
//    diag_CTC.setFromTriplets(triplet_list.begin(), triplet_list.end());
//    SpMat A2 = scaled_adjacency_mat_.transpose() * diag_CTC * scaled_adjacency_mat_;
//    for (int i = 0; i < num_iter_; ++i) {
//        std::cout << "Iteration " << i + 1 << std::endl;
////        std::cout << "Number o vertex " << num_vertex <<std::endl;
////        std::cout << "Dimension of matrix U = " << U.rows() << " " << U.cols() << std::endl;
////        std::cout << "Dimension of D = " << diagonal_vec.rows() << std::endl;
//        ProbabilSpectralDecomposition(U, diagonal_vec, scaled_adjacency_mat_, rank_, oversampling_p_, A1, A2);
//    }
//    VectorXd identity_vec(num_vertex);
//    for (int i = 0; i < num_vertex; ++i)
//        identity_vec(i) = 1;
//    SpMat scaled_adj_mat_product = scaled_adjacency_mat_.transpose() * scaled_adjacency_mat_;
//    simrank_ = MatrixXd(scaled_adj_mat_product);
//    simrank_.diagonal() = identity_vec;
//    simrank_ += U * diagonal_vec.asDiagonal() * U.transpose();
//    printf("Compute low rank approximation of SimRank... Done\n");
//}
//
//void SimRankLowrank::ProbabilSpectralDecomposition(MatrixXd& U, VectorXd& diag_D, const SpMat& scaled_adj_matrix, 
//                                                   int rank, int oversample_par, const SpMat& A1, 
//                                                   const SpMat& A2) {
//    // Generate matrix with elements from standard normal distribution
//    int num_vertex = scaled_adj_matrix.cols();
//    MatrixXd Z(num_vertex, rank + oversample_par);
//    std::random_device rd;
//    std::mt19937 generator(rd());
//    std::normal_distribution<> distribution(0.0, 1.0);
//    for (int i = 0; i < num_vertex; ++i) {
//        for (int j = 0; j < rank + oversample_par; ++j)
//            Z(i, j) = distribution(generator);
//    }
//    
////    SpMat CTC = scaled_adj_matrix.transpose() * scaled_adj_matrix;
////    SpMat A1 = scaled_adj_matrix.transpose() * CTC * scaled_adj_matrix;
////    MatrixXd D = MatrixXd(CTC).diagonal().asDiagonal();
////    MatrixXd A2 = scaled_adj_matrix.transpose() * D * scaled_adj_matrix;
////    MatrixXd A3 = scaled_adj_matrix.transpose() * U * diag_D.asDiagonal() * U.transpose() * scaled_adj_matrix;
////    MatrixXd A4 = (MatrixXd(A1) - A2 + A3).diagonal().asDiagonal();
//    // Form matrix Y = MZ
//    MatrixXd Y(num_vertex, rank + oversample_par);
////    MatrixXd M = MatrixXd(A1) - A2 + A3 - A4;
//    for (int i = 0; i < rank + oversample_par; ++i)
//        Y.col(i) = matvec(Z.col(i), U, diag_D, A1, A2);
//    
//    Eigen::JacobiSVD<MatrixXd> svd(Y, Eigen::ComputeThinU);
//    MatrixXd U_Y = svd.matrixU();
//    int rank_Y = rank + oversample_par; // FIXIT!
//    MatrixXd Q(num_vertex, rank_Y);
//    for (int i = 0; i < rank_Y; ++i)
//        Q.col(i) = U_Y.col(i);
//    
//    MatrixXd MQ(num_vertex, rank_Y);
//    for (int i = 0; i < rank_Y; ++i)
//        MQ.col(i) = matvec(Q.col(i), U, diag_D, A1, A2);
//    
//    MatrixXd B = Q.transpose() * MQ;
//    Eigen::JacobiSVD<MatrixXd> svd2(B, Eigen::ComputeThinU);
//    U = Q * svd2.matrixU();
//    diag_D = svd2.singularValues();
//}
//
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