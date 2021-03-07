// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix discreteStochModel(
        double nu, double nuA, 
        NumericVector probE,
        NumericVector probEP,
        NumericVector probA,
        NumericVector probP,
        NumericVector probI1,
        NumericVector probI1H,
        NumericVector probI1D,
        NumericVector probI2,
        NumericVector probH,
        NumericVector probHD,
        int tstop, arma::imat u, arma::mat C) {
    
    int i = 0, j = 0, k = 0;
      
    // set up auxiliary matrix for counts
    int nclasses = u.n_rows;
    int nages = u.n_cols;
    arma::imat u1(nclasses, nages);
    for(i = 0; i < nclasses; i++) {
        for(j = 0; j < nages; j++) {
            u1(i, j) = u(i, j);
        }
    }
    
    // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
    //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    
    // set up vector of number of infectives
    arma::vec uinf(nages);
    
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
    IntegerVector pathE(3);
    NumericVector mprobsE(3);
    IntegerVector pathI1(4);
    NumericVector mprobsI1(4);
    IntegerVector pathH(3);
    NumericVector mprobsH(3);
    tstart++;
    
    while(tstart < tstop) {
    
        // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
        //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
                
        for(j = 0; j < nages; j++) {
                                
            // H out
            mprobsH[0] = probH[j] * (1.0 - probHD[j]);
            mprobsH[1] = probH[j] * probHD[j];
            mprobsH[2] = 1.0 - probH[j];
            rmultinom(u1(9, j), mprobsH.begin(), 3, pathH.begin());
            u1(9, j) -= (pathH[0] + pathH[1]);
            u1(10, j) += pathH[0];
            u1(11, j) += pathH[1];
                        
            // I2RI
            i = R::rbinom(u1(7, j), probI2[j]);
            u1(7, j) -= i;
            u1(8, j) += i;
            
            // I1 out            
            mprobsI1[0] = probI1[j] * probI1H[j];
            mprobsI1[1] = probI1[j] * (1.0 - probI1D[j] - probI1H[j]);
            mprobsI1[2] = probI1[j] * probI1D[j];
            mprobsI1[3] = 1.0 - probI1[j];
            rmultinom(u1(5, j), mprobsI1.begin(), 4, pathI1.begin());
            u1(5, j) -= (pathI1[0] + pathI1[1] + pathI1[2]);
            u1(9, j) += pathI1[0];
            u1(7, j) += pathI1[1];
            u1(6, j) += pathI1[2];
                        
            // PI1
            i = R::rbinom(u1(4, j), probP[j]);
            u1(4, j) -= i;
            u1(5, j) += i;
                        
            // ARA
            i = R::rbinom(u1(2, j), probA[j]);
            u1(2, j) -= i;
            u1(3, j) += i;
                        
            // EA - EP
            mprobsE[0] = probE[j] * (1.0 - probEP[j]);
            mprobsE[1] = probE[j] * probEP[j];
            mprobsE[2] = 1.0 - probE[j];
            rmultinom(u1(1, j), mprobsE.begin(), 3, pathE.begin());
            u1(1, j) -= (pathE[0] + pathE[1]);
            u1(2, j) += pathE[0];
            u1(4, j) += pathE[1];
        }
        
        // update infective counts for rate
        for(i = 0; i < nages; i++) {
            uinf[i] = (double) nuA * u1(2, i) + nu * (u1(4, i) + u1(5, i) + u1(7, i));
        }
        
        // SE
        for(j = 0; j < nages; j++) {
            betaD = 0.7 * uinf / N;
            beta = C * betaD;
//            Rprintf("t = %d beta[%d] = %.25f\n", tstart, j, beta[j]);
            prob = 1 - exp(-beta(j, 0));
            i = R::rbinom(u1(0, j), prob);
            u1(0, j) -= i;
            u1(1, j) += i;
            
            betaN = 0.3 * uinf / N;
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

