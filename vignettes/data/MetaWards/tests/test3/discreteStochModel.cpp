// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix discreteStochModel(NumericVector pars, int tstop, arma::imat u, arma::mat C) {
    
    int i = 0, j = 0, k = 0;
  
    // extract parameters
    double nu = pars[0];
    double gammaE = pars[1];
    double gammaP = pars[2];
    double gammaI1 = pars[3];
      
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
        uinf[i] = (double) u(2, i) + u(3, i);
    }
    
    // extract population sizes
    arma::vec N (nages);
    for(i = 0; i < nages; i++) {
        N[i] = 0.0;
        for(j = 0; j < nclasses; j++) {
            N[i] += (double) u(j, i);
        }
    }
    
    // set up output
    int tstart = 0;
    IntegerMatrix out (tstop + 1, nclasses * nages + 1);
    out(0, 0) = 0;
    k = 1;
    for(i = 0; i < nclasses; i++) {
        for(j = 0; j < nages; j++) {
            out(0, k) = u(i, j);
            k++;
        }
    }
    
    // set up transmission rates
    arma::mat beta(nages, 1);
    beta = nu * C * uinf / N;
    
    // set up auxiliary variables
    double prob = 0.0;
    double probE = 1 - exp(-gammaE);
    double probP = 1 - exp(-gammaP);
    double probI1 = 1 - exp(-gammaI1);
    tstart++;
    
    while(tstart < tstop) {
        
        for(j = 0; j < nages; j++) {
            // SE
            prob = 1 - exp(-beta(j, 0));
            i = R::rbinom(u(0, j), prob);
            u1(0, j) -= i;
            u1(1, j) += i;
            
            // EP
            i = R::rbinom(u(1, j), probE);
            u1(1, j) -= i;
            u1(2, j) += i;
            
            // PI1
            i = R::rbinom(u(2, j), probP);
            u1(2, j) -= i;
            u1(3, j) += i;
            
            // I1DI
            i = R::rbinom(u(3, j), probI1);
            u1(3, j) -= i;
            u1(4, j) += i;
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
        
        // update counters
        for(i = 0; i < nclasses; i++) {
            for(j = 0; j < nages; j++) {
                u(i, j) = u1(i, j);
            }
        }
        for(i = 0; i < nages; i++) {
            uinf[i] = u(2, i) + u(3, i);
        }
        beta = nu * C * uinf / N;
        
        // update time 
        tstart++;
        Rprintf("tstart = %d\n", tstart);
    }
    Rprintf("hmm\n");
    return out;
}

