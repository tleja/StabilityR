#include <RcppArmadillo.h>
using namespace Rcpp;

// normalised variation of Information RcppArmadillo function.
// input: [x] partition matrix with n x m, where n is a number of nodes in the graph, and m is a number of partition vectors.
// output: square matrix m x m, with pair-wise variation of information between different partition vectors.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat varinfo(arma::mat x) {
	
	// number of nodes and a vector with node sequence 1:n
	int n = x.n_rows;
	arma::vec v = arma::ones(n);  
	arma::vec s = arma::linspace(1,n,n);
	
	// number of partition vectors, m
	int m = x.n_cols;
	arma::mat vimat = arma::zeros(m,m);
	
	for(int i=1; i < m; ++i) {
		
		// contigency matrix for partition vector, i
		arma::vec p1 = x.col(i);
		arma::mat ind1 = arma::join_cols(p1.t(), s.t());
		arma::umat indx1 = arma::conv_to<arma::umat>::from(ind1);
		arma::sp_mat a1(indx1-1, v);
		arma::sp_mat a11x = sum(a1,1);
		arma::mat a11 = arma::conv_to<arma::mat>::from(a11x);
	
		for(int j=0; j < i; ++j) {
	
			// contigency matrix for partition vector, j
			arma::vec p2 = x.col(j);
			arma::mat ind2 = arma::join_cols(s.t(),p2.t());
			arma::umat indx2 = arma::conv_to<arma::umat>::from(ind2);
			arma::sp_mat a2(indx2-1, v);
			arma::sp_mat a22x = sum(a2,0);
			arma::mat a22 = arma::conv_to<arma::mat>::from(a22x);

			// combined contingency matrix between i and j
			arma::sp_mat a12 = a1*a2;
			arma::mat aa = arma::conv_to<arma::mat>::from(a12);

			// vectors with values and indices exctracted from the combined contingency matrix
			arma::vec a_values = nonzeros(a12);
			arma::uvec a_indices = find(aa);
			arma::umat a_indmat = trans(ind2sub(size(aa), a_indices));
			arma::uvec a_row = a_indmat.col(0);
			arma::uvec a_col = a_indmat.col(1);
	
			// normalised variation of information between i and j
			arma::vec f1 = a11(a_row) % a22(a_col);
			arma::vec f2 = log(pow(a_values,2)/f1);
			double f3 = sum(a_values % f2);
			double vi = -1 / (n * log(n)) * f3;
			
			vimat(i,j) = vi;
		}
	}
	
	return vimat + trans(vimat);
}
