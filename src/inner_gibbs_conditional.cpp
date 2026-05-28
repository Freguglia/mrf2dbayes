#include <RcppArmadillo.h>
#include <numeric>   // std::iota
#include <cmath>     // std::exp
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
IntegerMatrix inner_gibbs_conditional(IntegerMatrix zinit, arma::fcube &cond_weights,
                                      IntegerMatrix R, arma::fcube &theta, int ncycles = 1){

    int N     = zinit.nrow();
    int M     = zinit.ncol();
    int latsize = N * M;
    int C     = (int)theta.n_rows - 1;
    int npos  = R.nrow();

    // Plain C++ arrays: no SEXP overhead, cache-friendly
    std::vector<double> lprobs(C + 1);
    std::vector<double> probs(C + 1);

    std::vector<int> order(latsize);
    std::iota(order.begin(), order.end(), 0);

    IntegerMatrix zout = Rcpp::clone(zinit);

    // Single RNG state sync for the entire function instead of
    // one pair per Rcpp::sample() call
    GetRNGstate();

    for(int t = 0; t < ncycles; t++){

        // Fisher-Yates shuffle via unif_rand()
        for(int i = latsize - 1; i > 0; i--){
            int j = (int)(unif_rand() * (i + 1));
            std::swap(order[i], order[j]);
        }

        for(int i = 0; i < latsize; i++){
            int x = order[i] / M;
            int y = order[i] % M;

            // Compute log-weights and track max inline
            double maxlp = -1e300;
            for(int k = 0; k <= C; k++){
                double H = (double)cond_weights(x, y, k);
                for(int r = 0; r < npos; r++){
                    int dx = R(r, 0), dy = R(r, 1);
                    int xp = x + dx, yp = y + dy;
                    int xm = x - dx, ym = y - dy;
                    if(xp >= 0 && xp < N && yp >= 0 && yp < M)
                        H += (double)theta(k, zout(xp, yp), r);
                    if(xm >= 0 && xm < N && ym >= 0 && ym < M)
                        H += (double)theta(zout(xm, ym), k, r);
                }
                lprobs[k] = H;
                if(H > maxlp) maxlp = H;
            }

            // Softmax and cumulative sum
            double total = 0.0;
            for(int k = 0; k <= C; k++){
                probs[k] = std::exp(lprobs[k] - maxlp);
                total += probs[k];
            }

            // Inline categorical sample — replaces Rcpp::sample(values, 1, ...)
            double u = unif_rand() * total;
            double cumsum = 0.0;
            int chosen = C;           // fallback to last category
            for(int k = 0; k < C; k++){
                cumsum += probs[k];
                if(u <= cumsum){ chosen = k; break; }
            }
            zout(x, y) = chosen;
        }
    }

    PutRNGstate();

    return zout;
}

