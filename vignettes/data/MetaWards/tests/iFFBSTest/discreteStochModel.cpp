// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix discreteStochModel(NumericVector pars, int tstop, arma::ivec u) {
    
    int i = 0, j = 0, k = 0;
  
    // extract parameters
    double nu = pars[0];
    double nuA = pars[1];
    double probE = pars[2];
    double probEP = pars[3];
    double probA = pars[4];
    double probP = pars[5];
    double probI1 = pars[6];
    double probI1H = pars[7];
    double probI1D = pars[8];
    double probI2 = pars[9];
    double probH = pars[10];
    double probHD = pars[11];
      
    // set up auxiliary matrix for counts
    int nclasses = u.n_elem;
    arma::ivec u1(nclasses);
    u1 = u;
    
    // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
    //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    
    // extract population sizes
    int N = 0;
    for(j = 0; j < nclasses; j++) {
        N += u1(j);
    }
    
    // set up output
    int tstart = 0;
    IntegerMatrix out(tstop, nclasses + 1);
    out(0, 0) = 0;
    for(i = 0; i < nclasses; i++) {
        out(0, i + 1) = u1(i);
    }
    
    // set up transmission rates
    double beta;
    
    // set up auxiliary variables
    double prob = 0.0;
    IntegerVector pathE(3);
    NumericVector mprobsE(3);
    mprobsE[0] = probE * (1.0 - probEP);
    mprobsE[1] = probE * probEP;
    mprobsE[2] = 1.0 - probE;
    IntegerVector pathI1(4);
    NumericVector mprobsI1(4);
    mprobsI1[0] = probI1 * probI1H;
    mprobsI1[1] = probI1 * (1.0 - probI1D - probI1H);
    mprobsI1[2] = probI1 * probI1D;
    mprobsI1[3] = 1.0 - probI1;
    IntegerVector pathH(3);
    NumericVector mprobsH(3);
    mprobsH[0] = probH * (1.0 - probHD);
    mprobsH[1] = probH * probHD;
    mprobsH[2] = 1.0 - probH;
    tstart++;
    
    while(tstart < tstop) {
    
        // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
        //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
                        
        // H out
        rmultinom(u1(9), mprobsH.begin(), 3, pathH.begin());
        u1(9) -= (pathH[0] + pathH[1]);
        u1(10) += pathH[0];
        u1(11) += pathH[1];
                    
        // I2RI
        i = R::rbinom(u1(7), probI2);
        u1(7) -= i;
        u1(8) += i;
        
        // I1 out
        rmultinom(u1(5), mprobsI1.begin(), 4, pathI1.begin());
        u1(5) -= (pathI1[0] + pathI1[1] + pathI1[2]);
        u1(9) += pathI1[0];
        u1(7) += pathI1[1];
        u1(6) += pathI1[2];
                    
        // PI1
        i = R::rbinom(u1(4), probP);
        u1(4) -= i;
        u1(5) += i;
                    
        // ARA
        i = R::rbinom(u1(2), probA);
        u1(2) -= i;
        u1(3) += i;
                    
        // EA - EP
        rmultinom(u1(1), mprobsE.begin(), 3, pathE.begin());
        u1(1) -= (pathE[0] + pathE[1]);
        u1(2) += pathE[0];
        u1(4) += pathE[1];
        
        prob = (double) nu * (nuA * u1(2) + u1(4) + u1(5) + u1(7));
        prob /= (double) N;
        prob = 1.0 - exp(-prob);
        
        i = R::rbinom(u1(0), prob);
        u1(0) -= i;
        u1(1) += i;
        
        // record output
        out(tstart, 0) = tstart;
        for(i = 0; i < nclasses; i++) {
            out(tstart, i + 1) = u1(i);
        }
        
        // update time 
        tstart++;
    }
    return out;
}

