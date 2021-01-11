#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::cube indicator_array(IntegerMatrix z, int C){

  int N = z.nrow();
	int M = z.ncol();
	arma::cube out = arma::zeros(N,M,C+1);

	for(int x = 0; x<N; x++){
		for(int y = 0; y<M; y++){
			out(x,y,z(x,y)) = 1;
		}
	}

	return(out);
}
