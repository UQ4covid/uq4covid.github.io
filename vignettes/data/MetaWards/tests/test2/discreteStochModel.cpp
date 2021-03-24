// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix discreteStochModel(NumericVector pars, int tstop, arma::imat u, arma::mat C) {
    
    int i = 0, j = 0, k = 0;
  
    // extract parameters
    double nu = pars[0];
    double probE = pars[1];
    double probP = pars[2];
    double probI1 = pars[3];
      
    // set up auxiliary matrix for counts
    int nclasses = u.n_rows;
    int nages = u.n_cols;
    arma::imat u1(nclasses, nages);
    for(i = 0; i < nclasses; i++) {
        for(j = 0; j < nages; j++) {
            u1(i, j) = u(i, j);
        }
    }
    
    // classes are: S, E, P, I1, DI
    
    // set up vector of number of infectives
    arma::vec uinf(nages);
    for(i = 0; i < nages; i++) {
        uinf[i] = (double) u1(2, i) + u1(3, i);
    }
    
    // extract population sizes
    arma::vec N (nages);
    for(i = 0; i < nages; i++) {
        N[i] = 0.0;
        for(j = 0; j < nclasses; j++) {
            N[i] += (double) u1(j, i);
        }
    }
    
    // set up output
    int tstart = 0;
    IntegerMatrix out (tstop + 1, nclasses * nages + 1);
    out(0, 0) = 0;
    k = 1;
    for(i = 0; i < nclasses; i++) {
        for(j = 0; j < nages; j++) {
            out(0, k) = u1(i, j);
            k++;
        }
    }
    
    // set up transmission rates
    arma::mat beta(nages, 1);
    arma::mat betaD(nages, 1);
    arma::mat betaN(nages, 1);
    
    // set up auxiliary variables
    double prob = 0.0;
    tstart++;
    
    while(tstart < tstop) {
        
        // update infective counts for rate
        for(i = 0; i < nages; i++) {
            uinf[i] = (double) u1(2, i) + u1(3, i);
        }
        
        for(j = 0; j < nages; j++) {
            
            // I1DI
            i = R::rbinom(u1(3, j), probI1);
            u1(3, j) -= i;
            u1(4, j) += i;
            
            // PI1
            i = R::rbinom(u1(2, j), probP);
            u1(2, j) -= i;
            u1(3, j) += i;
            
            // EP
            i = R::rbinom(u1(1, j), probE);
            u1(1, j) -= i;
            u1(2, j) += i;
        }
        
        // SE
        for(j = 0; j < nages; j++) {
            betaD = 0.7 * nu * uinf / N;
            beta = C * betaD;
//            Rprintf("t = %d beta[%d] = %.25f\n", tstart, j, beta[j]);
            prob = 1 - exp(-beta(j, 0));
            i = R::rbinom(u1(0, j), prob);
            u1(0, j) -= i;
            u1(1, j) += i;
            
            betaN = 0.3 * nu * uinf / N;
            beta = C * betaN;
//            Rprintf("t = %d beta[%d] = %.25f\n", tstart, j, beta[j]);
            prob = 1 - exp(-beta(j, 0));
            i = R::rbinom(u1(0, j), prob);
            u1(0, j) -= i;
            u1(1, j) += i;
        }
        
        // record output
        out(tstart, 0) = tstart;
        k = 1;
        for(i = 0; i < nclasses; i++) {
            for(j = 0; j < nages; j++) {
                out(tstart, k) = u1(i, j);
                k++;
            }
        }
        
        // update time 
        tstart++;
    }
    return out;
}

