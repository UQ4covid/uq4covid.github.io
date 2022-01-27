// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix discreteStochModel(NumericVector pars, int tstart, int tstop, NumericVector ageProbs,
                                 arma::imat workMoves, arma::mat C) {
    
    // workMoves is a matrix with columns: LADfrom, LADto, nmoves, nseeds
    
    // set up auxiliary matrix for counts
    int nclasses = 12;
    int nages = ageProbs.size();
    int nlads = (int) max(workMoves.col(0));
    int k;
    arma::uword i, j, l;
    
    // create cube object to store nclass x nage x nmove information
    arma::icube u1(nclasses, nages, workMoves.n_rows);
    arma::icube u1_day(nclasses, nages, nlads);
    arma::icube u1_night(nclasses, nages, nlads);
    u1.zeros();
    u1_day.zeros();
    u1_night.zeros();
    arma::imat N_day(nages, nlads);
    arma::imat N_night(nages, nlads);
    N_day.zeros();
    N_night.zeros();
    arma::mat pinf(nages, nlads);
    arma::imat origE(nages, workMoves.n_rows);
    
    // split workers and seed into age-classes probabilistically
    IntegerVector pathAge(nages);
    NumericVector seedProbs(nages);
    double nseeds = 0.0;
    for(i = 0; i < workMoves.n_rows; i++) {
        // split population by age
        rmultinom(workMoves(i, 2), ageProbs.begin(), nages, pathAge.begin());
        for(j = 0; j < nages; j++) {
            u1(0, j, i) = pathAge(j);
            // if(u1(0, j, i) > 1000000) Rprintf("u1 = %d i = %d j = %d\n", u1(0, j, i), i, j);
        }
        // set seeds within population (taking care to account for finite population size)
        for(l = 0; l < workMoves(i, 3); l++) {
            nseeds = 0.0;
            for(j = 0; j < nages; j++) {
                nseeds += u1(0, j, i);
            }
            for(j = 0; j < nages; j++) {
                seedProbs(j) = u1(0, j, i) / nseeds;
            }
            rmultinom(1, seedProbs.begin(), nages, pathAge.begin());
            for(j = 0; j < nages; j++) {
                u1(0, j, i) -= pathAge(j);
                u1(1, j, i) += pathAge(j);
                // if(u1(0, j, i) > 1000000) Rprintf("u1 = %d i = %d j = %d\n", u1(0, j, i), i, j);
                // if(u1(0, j, i) < 0) Rprintf("neg wk = %d u1 = %d i = %d j = %d\n", workMoves(i, 3), u1(0, j, i), i, j);
            }
        }
        // if(i == 196) {
        //     for(j = 0; j < nages; j++) {
        //         Rprintf("u1(0, %d, %d) = %d ", j, i, u1(0, j, i));
        //         Rprintf("u1(1, %d, %d) = %d\n", j, i, u1(1, j, i));
        //     }
        // }
    }
    // construct day and night counts
    for(l = 0; l < workMoves.n_rows; l++) {
        for(i = 0; i < nclasses; i++) {
            for(j = 0; j < nages; j++) {
                u1_day(i, j, (arma::uword) workMoves(l, 1) - 1) += u1(i, j, l);
                u1_night(i, j, (arma::uword) workMoves(l, 0) - 1) += u1(i, j, l);
            }
        }
        // if(l == 196) {
        //     for(j = 0; j < nages; j++) {
        //         Rprintf("u1_night(0, %d, %d) = %d ", j, workMoves(l, 0) - 1, u1_night(0, j, (arma::uword) workMoves(l, 0) - 1));
        //         Rprintf("u1_night(1, %d, %d) = %d\n", j, workMoves(l, 0) - 1, u1_night(1, j, (arma::uword) workMoves(l, 0) - 1));
        //         Rprintf("u1_day(0, %d, %d) = %d ", j, workMoves(l, 1) - 1, u1_day(0, j, (arma::uword) workMoves(l, 1) - 1));
        //         Rprintf("u1_day(1, %d, %d) = %d\n", j, workMoves(l, 1) - 1, u1_day(1, j, (arma::uword) workMoves(l, 1) - 1));
        //     }
        // }
    }
    for(i = 0; i < nlads; i++) {
        for(j = 0; j < nages; j++) {
            for(l = 0; l < nclasses; l++) {
                N_day(j, i) += u1_day(l, j, i);
                N_night(j, i) += u1_night(l, j, i);
            }
        }
    }
    for(j = 0; j < nages; j++) {
        if(sum(N_day.row(j)) != sum(N_night.row(j))) {
            Rprintf("Error in counts\n");
        }
    }
    
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
        for(i = 0; i < workMoves.n_rows; i++) {
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
                k = R::rbinom(u1(0, j, i), pinf(j, (arma::uword) workMoves(i, 1) - 1));
                // if(k > 0 && workMoves(i, 0) != 1 && workMoves(i, 1) != 1) Rprintf("hmmm\n");
                // if(k > 0) Rprintf("hmmm\n");
                // if(k > 100) Rprintf("day i = %d j = %d, u = %d k = %d\n", i, j, u1(0, j, i), k);
                // if(k < 0) Rprintf("kneg day i = %d j = %d, u = %d k = %d\n", i, j, u1(0, j, i), k);
                u1(0, j, i) -= k;
                u1(1, j, i) += k;
                u1_day(0, j, (arma::uword) workMoves(i, 1) - 1) -= k;
                u1_day(1, j, (arma::uword) workMoves(i, 1) - 1) += k;
                u1_night(0, j, (arma::uword) workMoves(i, 0) - 1) -= k;
                u1_night(1, j, (arma::uword) workMoves(i, 0) - 1) += k;
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
                k = R::rbinom(u1(0, j, i), pinf(j, (arma::uword) workMoves(i, 0) - 1));
                // if(k > 0 && workMoves(i, 0) != 1 && workMoves(i, 1) != 1) Rprintf("hmmm\n");
                // if(k > 0) Rprintf("hmmm\n");
                // if(k > 100) Rprintf("night i = %d j = %d, u = %d k = %d\n", i, j, u1(0, j, i), k);
                // if(k < 0) Rprintf("kneg night i = %d j = %d, u = %d k = %d\n", i, j, u1(0, j, i), k);
                u1(0, j, i) -= k;
                u1(1, j, i) += k;
                u1_day(0, j, (arma::uword) workMoves(i, 1) - 1) -= k;
                u1_day(1, j, (arma::uword) workMoves(i, 1) - 1) += k;
                u1_night(0, j, (arma::uword) workMoves(i, 0) - 1) -= k;
                u1_night(1, j, (arma::uword) workMoves(i, 0) - 1) += k;
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
            for(i = 0; i < workMoves.n_rows; i++) {
                
                // H out
                rmultinom(u1(9, j, i), mprobsH.begin(), 3, pathH.begin());
                u1(9, j, i) -= (pathH(0) + pathH(1));
                u1(10, j, i) += pathH(0);
                u1(11, j, i) += pathH(1);
                u1_day(9, j, (arma::uword) workMoves(i, 1) - 1) -= (pathH(0) + pathH(1));
                u1_day(10, j, (arma::uword) workMoves(i, 1) - 1) += pathH(0);
                u1_day(11, j, (arma::uword) workMoves(i, 1) - 1) += pathH(1);
                u1_night(9, j, (arma::uword) workMoves(i, 0) - 1) -= (pathH(0) + pathH(1));
                u1_night(10, j, (arma::uword) workMoves(i, 0) - 1) += pathH(0);
                u1_night(11, j, (arma::uword) workMoves(i, 0) - 1) += pathH(1);
                
                // I2RI
                k = R::rbinom(u1(7, j, i), probI2(j));
                u1(7, j, i) -= k;
                u1(8, j, i) += k;
                u1_day(7, j, (arma::uword) workMoves(i, 1) - 1) -= k;
                u1_day(8, j, (arma::uword) workMoves(i, 1) - 1) += k;
                u1_night(7, j, (arma::uword) workMoves(i, 0) - 1) -= k;
                u1_night(8, j, (arma::uword) workMoves(i, 0) - 1) += k;
                
                // I1 out
                rmultinom(u1(5, j, i), mprobsI1.begin(), 4, pathI1.begin());
                u1(5, j, i) -= (pathI1(0) + pathI1(1) + pathI1(2));
                u1(9, j, i) += pathI1(0);
                u1(7, j, i) += pathI1(1);
                u1(6, j, i) += pathI1(2);
                u1_day(5, j, (arma::uword) workMoves(i, 1) - 1) -= (pathI1(0) + pathI1(1) + pathI1(2));
                u1_day(9, j, (arma::uword) workMoves(i, 1) - 1) += pathI1(0);
                u1_day(7, j, (arma::uword) workMoves(i, 1) - 1) += pathI1(1);
                u1_day(6, j, (arma::uword) workMoves(i, 1) - 1) += pathI1(2);
                u1_night(5, j, (arma::uword) workMoves(i, 0) - 1) -= (pathI1(0) + pathI1(1) + pathI1(2));
                u1_night(9, j, (arma::uword) workMoves(i, 0) - 1) += pathI1(0);
                u1_night(7, j, (arma::uword) workMoves(i, 0) - 1) += pathI1(1);
                u1_night(6, j, (arma::uword) workMoves(i, 0) - 1) += pathI1(2);
                
                // PI1
                k = R::rbinom(u1(4, j, i), probP(j));
                u1(4, j, i) -= k;
                u1(5, j, i) += k;
                u1_day(4, j, (arma::uword) workMoves(i, 1) - 1) -= k;
                u1_day(5, j, (arma::uword) workMoves(i, 1) - 1) += k;
                u1_night(4, j, (arma::uword) workMoves(i, 0) - 1) -= k;
                u1_night(5, j, (arma::uword) workMoves(i, 0) - 1) += k;
                
                // ARA
                k = R::rbinom(u1(2, j, i), probA(j));
                u1(2, j, i) -= k;
                u1(3, j, i) += k;
                u1_day(2, j, (arma::uword) workMoves(i, 1) - 1) -= k;
                u1_day(3, j, (arma::uword) workMoves(i, 1) - 1) += k;
                u1_night(2, j, (arma::uword) workMoves(i, 0) - 1) -= k;
                u1_night(3, j, (arma::uword) workMoves(i, 0) - 1) += k;
                
                // E out
                rmultinom(origE(j, i), mprobsE.begin(), 3, pathE.begin());
                u1(1, j, i) -= (pathE(0) + pathE(1));
                u1(2, j, i) += pathE(0);
                u1(4, j, i) += pathE(1);
                u1_day(1, j, (arma::uword) workMoves(i, 1) - 1) -= (pathE(0) + pathE(1));
                u1_day(2, j, (arma::uword) workMoves(i, 1) - 1) += pathE(0);
                u1_day(4, j, (arma::uword) workMoves(i, 1) - 1) += pathE(1);
                u1_night(1, j, (arma::uword) workMoves(i, 0) - 1) -= (pathE(0) + pathE(1));
                u1_night(2, j, (arma::uword) workMoves(i, 0) - 1) += pathE(0);
                u1_night(4, j, (arma::uword) workMoves(i, 0) - 1) += pathE(1);
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

