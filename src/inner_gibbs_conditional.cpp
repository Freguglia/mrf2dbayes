#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix inner_gibbs_conditional(IntegerMatrix zinit, arma::fcube &cond_weights,
                                      IntegerMatrix R, arma::fcube &theta, int ncycles = 1){

    int N = zinit.nrow();
    int M = zinit.ncol();
    int latsize = N*M;
    int C = theta.n_rows - 1;
    int x,y;
    int dx,dy;
    int npos = R.nrow();
    double H;
  
    NumericVector lprobs(C+1);
    NumericVector probs(C+1);
    IntegerVector values = seq_len(C+1) - 1;
    
    IntegerVector order = seq_len(latsize) - 1;
    IntegerMatrix zout = Rcpp::clone(zinit);

    for(int t=0; t<ncycles; t++){
        order = sample(order, latsize, false);
        for(int i=0; i<latsize; i++){
            x = order[i] / M; y = order[i] % M;
            for(int k=0; k<=C; k++){
                H = 0.0;
                for(int r=0; r<npos; r++){
                    dx = R(r,0); dy = R(r,1);
                    if(x + dx < N && x + dx >= 0 && y + dy < M && y + dy >= 0){
                        H = H + theta(k, zout(x+dx, y+dy), r);
                    }
                    if(x - dx < N && x - dx >= 0 && y - dy < M && y - dy >= 0){
                        H = H + theta(zout(x-dx,y-dy), k, r);
                    }
                }
                lprobs[k] = H + cond_weights(x,y,k);
            }
            lprobs = lprobs - max(lprobs);
            probs = exp(lprobs);
            zout(x,y) = sample(values, 1, false, probs)[0];
        }
    }

    return(zout);
}

