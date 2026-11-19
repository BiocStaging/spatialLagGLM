#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//' Compute Row-Normalized Spatial Weights
//' @param coords Numeric matrix of physical spots coordinates.
//' @param radius Distance radius.
//' @return Matrix of physical distance weights.
//' @keywords internal
//' @noRd
// [[Rcpp::export]]
arma::mat rcpp_spatial_weights(arma::mat coords, double radius) {
    int n = coords.n_rows;
    arma::mat W = arma::zeros<arma::mat>(n, n);
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dist = arma::norm(coords.row(i) - coords.row(j), 2);
            if (dist <= radius && dist > 0) {
                W(i, j) = 1.0 / dist;
                W(j, i) = 1.0 / dist;
            }
        }
    }
    
    for (int i = 0; i < n; i++) {
        double row_sum = arma::sum(W.row(i));
        if (row_sum > 0) {
            W.row(i) /= row_sum;
        }
    }
    
    return W;
}
