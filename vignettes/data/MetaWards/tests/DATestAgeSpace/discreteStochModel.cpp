// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix discreteStochModel(NumericVector pars, int tstart, int tstop, 
                                 arma::imat u1_moves, arma::icube u1, 
                                 arma::icube u1_day, arma::icube u1_night,
                                 arma::imat N_day, arma::imat N_night,
                                 arma::mat C) {
    
    // u1_moves is a matrix with columns: LADfrom, LADto
    // u1 is a 3D array with dimensions: nclasses x nages x nmoves
    //          each row of u1 must match u1_moves
    // u1_day and u1_night are 3D arrays with dimensions: nclasses x nages x nlads
    // N_day and N_night are nages x nlad matrices of population counts
    
    // set up auxiliary matrix for counts
    int nclasses = u1_day.n_rows;
    int nages = u1_day.n_cols;
    int nlads = u1_day.n_slices;
    int k;
    arma::uword i, j, l;
    
    // create objects to store move information
    arma::mat pinf(nages, nlads);
    arma::imat origE(nages, u1.n_slices);
    
    // extract parameters
    double nu = pars(0);
    double nuA = pars(1);
    arma::vec probE(nages);
    arma::vec probEP(nages);
    arma::vec probA(nages);
    arma::vec probP(nages);
    arma::vec probI1(nages);
    arma::vec probI1H(nages);
    arma::vec probI1D(nages);
    arma::vec probI2(nages);
    arma::vec probH(nages);
    arma::vec probHD(nages);
    
    for(j = 0; j < nages; j++) {
        probE(j) = pars(j + 2);
        probEP(j) = pars(j + nages + 2);
        probA(j) = pars(j + 2 * nages + 2);
        probP(j) = pars(j + 3 * nages + 2);
        probI1(j) = pars(j + 4 * nages + 2);
        probI1H(j) = pars(j + 5 * nages + 2);
        probI1D(j) = pars(j + 6 * nages + 2);
        probI2(j) = pars(j + 7 * nages + 2);
        probH(j) = pars(j + 8 * nages + 2);
        probHD(j) = pars(j + 9 * nages + 2);
    }
    
    // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
    //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    
    // // set up output
    // Rprintf("new\n");
    // for(j = 0; j < nages; j++) {
    //     Rprintf("u1(0, %d, 0) = %d ", j, u1_night(0, j, 0));
    //     Rprintf("u1(1, %d, 0) = %d\n", j, u1_night(1, j, 0));
    // }
    // Rprintf("new\n");
    // for(j = 0; j < nages; j++) {
    //     Rprintf("u1(0, %d, 282) = %d ", j, u1_day(0, j, 282));
    //     Rprintf("u1(1, %d, 282) = %d\n", j, u1_day(1, j, 282));
    // }
    
    IntegerMatrix out(tstop - tstart + 1, nclasses * nages * nlads + 1);
    int tcurr = 0;
    out(0, 0) = tstart;
    k = 1;
    for(i = 0; i < nclasses; i++) {
        for(j = 0; j < nages; j++) {
            for(l = 0; l < nlads; l++) {
                out(0, k) = u1_night(i, j, l);
                k++;
            }
        }
    }
    
    // set up vector of number of infectives
    arma::vec uinf(nages);
    
    // set up transmission rates
    arma::mat beta(nages, 1);
    
    // set up auxiliary variables
    double prob = 0.0;
    IntegerVector pathE(3);
    NumericVector mprobsE(3);
    IntegerVector pathI1(4);
    NumericVector mprobsI1(4);
    IntegerVector pathH(3);
    NumericVector mprobsH(3);
    tcurr++;
    tstart++;
    
    while(tstart <= tstop) {
    
        // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
        //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
        
        // save current number of infectives for later transitions
        for(i = 0; i < u1.n_slices; i++) {
            for(j = 0; j < nages; j++) {
                origE(j, i) = u1(1, j, i);
            }
        }
        
        // transmission probabilities (day), loop over LADs
        for(i = 0; i < u1_day.n_slices; i++) {
        
            // update infective counts for rate
            for(j = 0; j < nages; j++) {
                uinf(j) = (double) nuA * u1_day(2, j, i) + nu * (u1_day(4, j, i) + u1_day(5, j, i) + u1_day(7, j, i));
            }
        
            // SE
            beta = 0.7 * C * (uinf / N_day.col(i));
            for(j = 0; j < nages; j++) {
                // if(beta(j, 0) > 0.0) Rprintf("t = %d beta[%d] = %.25f\n", tstart, j, beta(j, 0));
                pinf(j, i) = 1.0 - exp(-beta(j, 0));
            }
        }
        // transmission events (day), loop over network
        for(i = 0; i < u1.n_slices; i++) {
            for(j = 0; j < nages; j++) {
                // if(u1(0, j, i) > 1000000) Rprintf("t = %d u1 = %d i = %d j = %d\n", tstart, u1(0, j, i), i, j);
                k = R::rbinom(u1(0, j, i), pinf(j, (arma::uword) u1_moves(i, 1) - 1));
                // if(k > 0 && u1_moves(i, 0) != 1 && u1_moves(i, 1) != 1) Rprintf("hmmm\n");
                // if(k > 0) Rprintf("hmmm\n");
                // if(k > 100) Rprintf("day i = %d j = %d, u = %d k = %d\n", i, j, u1(0, j, i), k);
                // if(k < 0) Rprintf("kneg day i = %d j = %d, u = %d k = %d\n", i, j, u1(0, j, i), k);
                u1(0, j, i) -= k;
                u1(1, j, i) += k;
                u1_day(0, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(1, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(0, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(1, j, (arma::uword) u1_moves(i, 0) - 1) += k;
            }
        }
        
        // transmission probabilities (night), loop over LADs
        for(i = 0; i < u1_night.n_slices; i++) {
            // Rprintf("N_Night(%d) = %d\n", i, N_night(i));
            
            // update infective counts for rate
            for(j = 0; j < nages; j++) {
                uinf(j) = (double) nuA * u1_night(2, j, i) + nu * (u1_night(4, j, i) + u1_night(5, j, i) + u1_night(7, j, i));
                // if(uinf(j) > 0) Rprintf("uinf(j) = %d\n", uinf(j));
            }
            
            // SE
            beta = 0.3 * C * (uinf / N_night.col(i));
            for(j = 0; j < nages; j++) {
                // if(beta(j, 0) > 0.0) Rprintf("night t = %d beta[%d] = %.25f\n", tstart, j, beta(j, 0));
                pinf(j, i) = 1.0 - exp(-beta(j, 0));
            }
        }
        // transmission events (night), loop over network
        for(i = 0; i < u1.n_slices; i++) {
            for(j = 0; j < nages; j++) {
                // if(u1(0, j, i) > 1000000) Rprintf("t = %d u1 = %d i = %d j = %d\n", tstart, u1(0, j, i), i, j);
                k = R::rbinom(u1(0, j, i), pinf(j, (arma::uword) u1_moves(i, 0) - 1));
                // if(k > 0 && u1_moves(i, 0) != 1 && u1_moves(i, 1) != 1) Rprintf("hmmm\n");
                // if(k > 0) Rprintf("hmmm\n");
                // if(k > 100) Rprintf("night i = %d j = %d, u = %d k = %d\n", i, j, u1(0, j, i), k);
                // if(k < 0) Rprintf("kneg night i = %d j = %d, u = %d k = %d\n", i, j, u1(0, j, i), k);
                u1(0, j, i) -= k;
                u1(1, j, i) += k;
                u1_day(0, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(1, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(0, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(1, j, (arma::uword) u1_moves(i, 0) - 1) += k;
            }
        }
            
        // conduct transition moves, loop over age classes                
        for(j = 0; j < nages; j++) {
            
            // (These probs could be pre-calculated an stored as matrix as long
            // as rmultinom can use rows of matrices as inputs?)
            
            // transition probs out of hospital
            mprobsH(0) = probH(j) * (1.0 - probHD(j));
            mprobsH(1) = probH(j) * probHD(j);
            mprobsH(2) = 1.0 - probH(j);
            
            // transition probs out of I1
            mprobsI1(0) = probI1(j) * probI1H(j);
            mprobsI1(1) = probI1(j) * (1.0 - probI1D(j) - probI1H(j));
            mprobsI1(2) = probI1(j) * probI1D(j);
            mprobsI1(3) = 1.0 - probI1(j);
            
            // transition probs out of E
            mprobsE(0) = probE(j) * (1.0 - probEP(j));
            mprobsE(1) = probE(j) * probEP(j);
            mprobsE(2) = 1.0 - probE(j);
            
            // loop over network
            for(i = 0; i < u1.n_slices; i++) {
                
                // H out
                rmultinom(u1(9, j, i), mprobsH.begin(), 3, pathH.begin());
                u1(9, j, i) -= (pathH(0) + pathH(1));
                u1(10, j, i) += pathH(0);
                u1(11, j, i) += pathH(1);
                u1_day(9, j, (arma::uword) u1_moves(i, 1) - 1) -= (pathH(0) + pathH(1));
                u1_day(10, j, (arma::uword) u1_moves(i, 1) - 1) += pathH(0);
                u1_day(11, j, (arma::uword) u1_moves(i, 1) - 1) += pathH(1);
                u1_night(9, j, (arma::uword) u1_moves(i, 0) - 1) -= (pathH(0) + pathH(1));
                u1_night(10, j, (arma::uword) u1_moves(i, 0) - 1) += pathH(0);
                u1_night(11, j, (arma::uword) u1_moves(i, 0) - 1) += pathH(1);
                
                // I2RI
                k = R::rbinom(u1(7, j, i), probI2(j));
                u1(7, j, i) -= k;
                u1(8, j, i) += k;
                u1_day(7, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(8, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(7, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(8, j, (arma::uword) u1_moves(i, 0) - 1) += k;
                
                // I1 out
                rmultinom(u1(5, j, i), mprobsI1.begin(), 4, pathI1.begin());
                u1(5, j, i) -= (pathI1(0) + pathI1(1) + pathI1(2));
                u1(9, j, i) += pathI1(0);
                u1(7, j, i) += pathI1(1);
                u1(6, j, i) += pathI1(2);
                u1_day(5, j, (arma::uword) u1_moves(i, 1) - 1) -= (pathI1(0) + pathI1(1) + pathI1(2));
                u1_day(9, j, (arma::uword) u1_moves(i, 1) - 1) += pathI1(0);
                u1_day(7, j, (arma::uword) u1_moves(i, 1) - 1) += pathI1(1);
                u1_day(6, j, (arma::uword) u1_moves(i, 1) - 1) += pathI1(2);
                u1_night(5, j, (arma::uword) u1_moves(i, 0) - 1) -= (pathI1(0) + pathI1(1) + pathI1(2));
                u1_night(9, j, (arma::uword) u1_moves(i, 0) - 1) += pathI1(0);
                u1_night(7, j, (arma::uword) u1_moves(i, 0) - 1) += pathI1(1);
                u1_night(6, j, (arma::uword) u1_moves(i, 0) - 1) += pathI1(2);
                
                // PI1
                k = R::rbinom(u1(4, j, i), probP(j));
                u1(4, j, i) -= k;
                u1(5, j, i) += k;
                u1_day(4, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(5, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(4, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(5, j, (arma::uword) u1_moves(i, 0) - 1) += k;
                
                // ARA
                k = R::rbinom(u1(2, j, i), probA(j));
                u1(2, j, i) -= k;
                u1(3, j, i) += k;
                u1_day(2, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(3, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(2, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(3, j, (arma::uword) u1_moves(i, 0) - 1) += k;
                
                // E out
                rmultinom(origE(j, i), mprobsE.begin(), 3, pathE.begin());
                u1(1, j, i) -= (pathE(0) + pathE(1));
                u1(2, j, i) += pathE(0);
                u1(4, j, i) += pathE(1);
                u1_day(1, j, (arma::uword) u1_moves(i, 1) - 1) -= (pathE(0) + pathE(1));
                u1_day(2, j, (arma::uword) u1_moves(i, 1) - 1) += pathE(0);
                u1_day(4, j, (arma::uword) u1_moves(i, 1) - 1) += pathE(1);
                u1_night(1, j, (arma::uword) u1_moves(i, 0) - 1) -= (pathE(0) + pathE(1));
                u1_night(2, j, (arma::uword) u1_moves(i, 0) - 1) += pathE(0);
                u1_night(4, j, (arma::uword) u1_moves(i, 0) - 1) += pathE(1);
            }
        }
        
        // record output
        out(tcurr, 0) = tstart;
        k = 1;
        for(i = 0; i < nclasses; i++) {
            for(j = 0; j < nages; j++) {
                for(l = 0; l < nlads; l++) {
                    out(tcurr, k) = u1_night(i, j, l);
                    k++;
                }
            }
        }
        
        // update time 
        tcurr++;
        tstart++;
    }
    return out;
}

