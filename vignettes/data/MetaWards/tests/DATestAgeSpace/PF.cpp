// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;

// log-sum-exp function to prevent numerical overflow
double log_sum_exp(arma::vec x, int mn = 0) {
    double maxx = max(x);
    double y = maxx + log(sum(exp(x - maxx)));
    if(mn == 1) y = y - log(x.n_elem);
    return y;
}

// log truncated Skellam density (note, no checks on inputs)
double ldtskellam_cpp(int x, double lambda1, double lambda2, int LB = 0, int UB = -1) {
    
    // if UB < LB then returns untruncated density
    
    // from "skellam" source code: log(besselI(y, nu)) == y + log(besselI(y, nu, TRUE))
    double ldens = -(lambda1 + lambda2) + (x / 2.0) * (log(lambda1) - log(lambda2)) +
        log(R::bessel_i(2.0 * sqrt(lambda1 * lambda2), abs(x), true)) + 2.0 * sqrt(lambda1 * lambda2);
    
    // if truncated then adjust log-density
    if(UB >= LB) {
        // declare variables
        int normsize = UB - LB + 1;
        if(normsize > 1) {
            // declare variables
            arma::vec norm (normsize);
            norm.zeros();
            double mnorm = 0.0, snorm = 0.0;
            for(arma::uword j = LB; j <= UB; j++) {
                norm(j - LB) = -(lambda1 + lambda2) + (j / 2.0) * (log(lambda1) - log(lambda2)) +
                    log(R::bessel_i(2.0 * sqrt(lambda1 * lambda2), abs(j), true)) + 2.0 * sqrt(lambda1 * lambda2);
            }
            ldens -= log_sum_exp(norm);
        } else {
            ldens = 0.0;
        }
    }
    return ldens;
}

// truncated Skellam sampler
int rtskellam_cpp(double lambda1, double lambda2, int LB = 0, int UB = -1) {
    
    // if UB < LB then returns untruncated density
    
    int x = 0;
    // if bounds are the same, then return the
    // only viable value
    if(LB == UB) {
        x = LB;
    } else {
        // draw untruncated sample
        x = R::rpois(lambda1) - R::rpois(lambda2);
        
        // rejection sample if necessary
        if(UB > LB) {
            int k = 0;
            int ntries = 1000;
            while((x < LB || x > UB) && k < ntries) {
                x = R::rpois(lambda1) - R::rpois(lambda2);
                k++;
            }
            // if rejection sampling doesn't work
            // then try long-hand way
            if(k == ntries) {
                // Rprintf("LB = %d UB = %d\n", LB, UB);
                double u = R::runif(0.0, 1.0);
                k = 0;
                double xdens = exp(ldtskellam_cpp(LB + k, lambda1, lambda2, LB, UB));
                while(u > xdens && k < (UB - LB + 1)) {
                    k++;
                    xdens += exp(ldtskellam_cpp(LB + k, lambda1, lambda2, LB, UB));
                }
                if(k == (UB - LB + 1) && u > xdens) {
                    stop("Something wrong in truncated Skellam sampling\n");
                }
                x = LB + k;
            }
        }
    }
    return x;
}

// simulation model
void discreteStochModel(arma::vec pars, int tstart, int tstop, 
                        arma::imat *u1_moves, arma::icube *u1, arma::icube *u1_day, arma::icube *u1_night,
                        arma::imat *N_day, arma::imat *N_night, arma::mat *pinf, arma::imat *origE, arma::mat C) {
    
    // u1_moves is a matrix with columns: LADfrom, LADto
    // u1 is a 3D array with dimensions: nclasses x nages x nmoves
    //          each row of u1 must match u1_moves
    // u1_day/night are 3D arrays with dimensions: nclasses x nages x nlads
    // N_day/night are nlad x nage matrices of population counts
    // pinf is nages x nlads auxiliary matrix
    // origE is nages x nmoves auxiliary matrix
    
    // set up auxiliary matrix for counts
    int nclasses = (*u1).n_rows;
    int nages = (*u1).n_cols;
    int nlads = max((*u1_moves).col(0));
    int k;
    arma::uword i, j;
    
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
    
    int tcurr = 0;
    
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
        for(i = 0; i < (*u1).n_slices; i++) {
            for(j = 0; j < nages; j++) {
                (*origE)(j, i) = (*u1)(1, j, i);
            }
        }
        
        // transmission probabilities (day), loop over LADs
        for(i = 0; i < (*u1_day).n_slices; i++) {
            
            // update infective counts for rate
            for(j = 0; j < nages; j++) {
                uinf(j) = (double) nuA * (*u1_day)(2, j, i) + nu * ((*u1_day)(4, j, i) + (*u1_day)(5, j, i) + (*u1_day)(7, j, i));
            }
            
            // SE
            beta = 0.7 * C * (uinf / (*N_day).col(i));
            for(j = 0; j < nages; j++) {
                (*pinf)(j, i) = 1.0 - exp(-beta(j, 0));
            }
        }
        // transmission events (day), loop over network
        for(i = 0; i < (*u1).n_slices; i++) {
            for(j = 0; j < nages; j++) {
                k = R::rbinom((*u1)(0, j, i), (*pinf)(j, (arma::uword) (*u1_moves)(i, 1) - 1));
                (*u1)(0, j, i) -= k;
                (*u1)(1, j, i) += k;
                (*u1_day)(0, j, (arma::uword) (*u1_moves)(i, 1) - 1) -= k;
                (*u1_day)(1, j, (arma::uword) (*u1_moves)(i, 1) - 1) += k;
                (*u1_night)(0, j, (arma::uword) (*u1_moves)(i, 0) - 1) -= k;
                (*u1_night)(1, j, (arma::uword) (*u1_moves)(i, 0) - 1) += k;
            }
        }
        
        // transmission probabilities (night), loop over LADs
        for(i = 0; i < (*u1_night).n_slices; i++) {
            
            // update infective counts for rate
            for(j = 0; j < nages; j++) {
                uinf(j) = (double) nuA * (*u1_night)(2, j, i) + nu * ((*u1_night)(4, j, i) + (*u1_night)(5, j, i) + (*u1_night)(7, j, i));
            }
            
            // SE
            beta = 0.3 * C * (uinf / (*N_night).col(i));
            for(j = 0; j < nages; j++) {
                (*pinf)(j, i) = 1.0 - exp(-beta(j, 0));
            }
        }
        // transmission events (night), loop over network
        for(i = 0; i < (*u1).n_slices; i++) {
            for(j = 0; j < nages; j++) {
                k = R::rbinom((*u1)(0, j, i), (*pinf)(j, (arma::uword) (*u1_moves)(i, 0) - 1));
                (*u1)(0, j, i) -= k;
                (*u1)(1, j, i) += k;
                (*u1_day)(0, j, (arma::uword) (*u1_moves)(i, 1) - 1) -= k;
                (*u1_day)(1, j, (arma::uword) (*u1_moves)(i, 1) - 1) += k;
                (*u1_night)(0, j, (arma::uword) (*u1_moves)(i, 0) - 1) -= k;
                (*u1_night)(1, j, (arma::uword) (*u1_moves)(i, 0) - 1) += k;
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
            for(i = 0; i < (*u1).n_slices; i++) {
                
                // H out
                rmultinom((*u1)(9, j, i), mprobsH.begin(), 3, pathH.begin());
                (*u1)(9, j, i) -= (pathH(0) + pathH(1));
                (*u1)(10, j, i) += pathH(0);
                (*u1)(11, j, i) += pathH(1);
                (*u1_day)(9, j, (arma::uword) (*u1_moves)(i, 1) - 1) -= (pathH(0) + pathH(1));
                (*u1_day)(10, j, (arma::uword) (*u1_moves)(i, 1) - 1) += pathH(0);
                (*u1_day)(11, j, (arma::uword) (*u1_moves)(i, 1) - 1) += pathH(1);
                (*u1_night)(9, j, (arma::uword) (*u1_moves)(i, 0) - 1) -= (pathH(0) + pathH(1));
                (*u1_night)(10, j, (arma::uword) (*u1_moves)(i, 0) - 1) += pathH(0);
                (*u1_night)(11, j, (arma::uword) (*u1_moves)(i, 0) - 1) += pathH(1);
                
                // I2RI
                k = R::rbinom((*u1)(7, j, i), probI2(j));
                (*u1)(7, j, i) -= k;
                (*u1)(8, j, i) += k;
                (*u1_day)(7, j, (arma::uword) (*u1_moves)(i, 1) - 1) -= k;
                (*u1_day)(8, j, (arma::uword) (*u1_moves)(i, 1) - 1) += k;
                (*u1_night)(7, j, (arma::uword) (*u1_moves)(i, 0) - 1) -= k;
                (*u1_night)(8, j, (arma::uword) (*u1_moves)(i, 0) - 1) += k;
                
                // I1 out
                rmultinom((*u1)(5, j, i), mprobsI1.begin(), 4, pathI1.begin());
                (*u1)(5, j, i) -= (pathI1(0) + pathI1(1) + pathI1(2));
                (*u1)(9, j, i) += pathI1(0);
                (*u1)(7, j, i) += pathI1(1);
                (*u1)(6, j, i) += pathI1(2);
                (*u1_day)(5, j, (arma::uword) (*u1_moves)(i, 1) - 1) -= (pathI1(0) + pathI1(1) + pathI1(2));
                (*u1_day)(9, j, (arma::uword) (*u1_moves)(i, 1) - 1) += pathI1(0);
                (*u1_day)(7, j, (arma::uword) (*u1_moves)(i, 1) - 1) += pathI1(1);
                (*u1_day)(6, j, (arma::uword) (*u1_moves)(i, 1) - 1) += pathI1(2);
                (*u1_night)(5, j, (arma::uword) (*u1_moves)(i, 0) - 1) -= (pathI1(0) + pathI1(1) + pathI1(2));
                (*u1_night)(9, j, (arma::uword) (*u1_moves)(i, 0) - 1) += pathI1(0);
                (*u1_night)(7, j, (arma::uword) (*u1_moves)(i, 0) - 1) += pathI1(1);
                (*u1_night)(6, j, (arma::uword) (*u1_moves)(i, 0) - 1) += pathI1(2);
                
                // PI1
                k = R::rbinom((*u1)(4, j, i), probP(j));
                (*u1)(4, j, i) -= k;
                (*u1)(5, j, i) += k;
                (*u1_day)(4, j, (arma::uword) (*u1_moves)(i, 1) - 1) -= k;
                (*u1_day)(5, j, (arma::uword) (*u1_moves)(i, 1) - 1) += k;
                (*u1_night)(4, j, (arma::uword) (*u1_moves)(i, 0) - 1) -= k;
                (*u1_night)(5, j, (arma::uword) (*u1_moves)(i, 0) - 1) += k;
                
                // ARA
                k = R::rbinom((*u1)(2, j, i), probA(j));
                (*u1)(2, j, i) -= k;
                (*u1)(3, j, i) += k;
                (*u1_day)(2, j, (arma::uword) (*u1_moves)(i, 1) - 1) -= k;
                (*u1_day)(3, j, (arma::uword) (*u1_moves)(i, 1) - 1) += k;
                (*u1_night)(2, j, (arma::uword) (*u1_moves)(i, 0) - 1) -= k;
                (*u1_night)(3, j, (arma::uword) (*u1_moves)(i, 0) - 1) += k;
                
                // E out
                rmultinom((*origE)(j, i), mprobsE.begin(), 3, pathE.begin());
                (*u1)(1, j, i) -= (pathE(0) + pathE(1));
                (*u1)(2, j, i) += pathE(0);
                (*u1)(4, j, i) += pathE(1);
                (*u1_day)(1, j, (arma::uword) (*u1_moves)(i, 1) - 1) -= (pathE(0) + pathE(1));
                (*u1_day)(2, j, (arma::uword) (*u1_moves)(i, 1) - 1) += pathE(0);
                (*u1_day)(4, j, (arma::uword) (*u1_moves)(i, 1) - 1) += pathE(1);
                (*u1_night)(1, j, (arma::uword) (*u1_moves)(i, 0) - 1) -= (pathE(0) + pathE(1));
                (*u1_night)(2, j, (arma::uword) (*u1_moves)(i, 0) - 1) += pathE(0);
                (*u1_night)(4, j, (arma::uword) (*u1_moves)(i, 0) - 1) += pathE(1);
            }
        }
        
        // update time 
        tcurr++;
        tstart++;
    }
    return;
}

// void test(arma::icube *u1) {
//     (*u1)(0, 0, 0) = 222;
//     return;
// }

// [[Rcpp::export]]
List PF_cpp (arma::vec pars, arma::mat C, arma::imat data, int nclasses, int nages, int nlads, arma::imat u1_moves, 
         arma::icube u1_comb, int ndays, int npart, int MD, double a1, double a2, double b, double a_dis, 
         double b_dis, int saveAll) {
    
    // set counters
    arma::uword i, j, l;
    int tempLB = 0;
    int t = 0;
    
    // split u1 up into different LADs
    std::vector<arma::icube> u1(npart);
    std::vector<arma::icube> u1_new(npart);
    
    // for(i = 0; i < 2; i++) {
    //     u1[i] = u1_comb;
    //     Rprintf("u1[%d](0, 0, 0) = %d\n", i, u1[i](0, 0, 0));
    // }
    // test(&u1[1]);
    // for(i = 0; i < 2; i++) {
    //     Rprintf("u1[%d](0, 0, 0) = %d\n", i, u1[i](0, 0, 0));
    // }
    // stop("STOP RIGHT NOW\n");
    
    // std::vector<arma::icube> cu(npart);
    std::vector<arma::icube> u1_night(npart);
    std::vector<arma::icube> u1_day(npart);
    std::vector<arma::icube> u1_night_new(npart);
    std::vector<arma::icube> u1_day_new(npart);
    arma::imat N_day(nages, nlads);
    arma::imat N_night(nages, nlads);
    arma::icube tu1_day(nclasses, nages, nlads);
    arma::icube tu1_night(nclasses, nages, nlads);
    tu1_day.zeros();
    tu1_night.zeros();
    N_day.zeros();
    N_night.zeros();
    for(i = 0; i < u1_moves.n_rows; i++) {
        for(j = 0; j < nages; j++) {
            for(l = 0; l < nclasses; l++) {
                tu1_day(l, j, (arma::uword) u1_moves(i, 1) - 1) += u1_comb(l, j, i);
                tu1_night(l, j, (arma::uword) u1_moves(i, 0) - 1) += u1_comb(l, j, i);
                N_day(j, (arma::uword) u1_moves(i, 1) - 1) += u1_comb(l, j, i);
                N_night(j, (arma::uword) u1_moves(i, 0) - 1) += u1_comb(l, j, i);
            }
        }
    }
    for(i = 0; i < npart; i++) {
        u1[i] = u1_comb;
        u1_new[i] = u1_comb;
        // cu[i] = u1_comb;
        u1_day[i] = tu1_day;
        u1_night[i] = tu1_night;
        u1_day_new[i] = tu1_day;
        u1_night_new[i] = tu1_night;
    }
    // auxiliary objects
    arma::mat pinf(nages, nlads);
    arma::imat origE(nages, u1_comb.n_slices);
    
    // set up weight vector
    NumericVector weights (npart);
    IntegerVector inds(npart);
    double wnorm = 0.0;
    
    // set up temporary objects
    int incsize = (MD == 1 ? u1_moves.n_rows:nlads);
    arma::imat DHinc (nages, incsize); DHinc.zeros();
    arma::imat RHinc (nages, incsize); RHinc.zeros();
    arma::imat Hinc (nages, incsize); Hinc.zeros();
    arma::imat DIinc (nages, incsize); DIinc.zeros();
    arma::imat RIinc (nages, incsize); RIinc.zeros();
    arma::imat I2inc (nages, incsize); I2inc.zeros();
    arma::imat I1inc (nages, incsize); I1inc.zeros();
    arma::imat Pinc (nages, incsize); Pinc.zeros();
    arma::imat RAinc (nages, incsize); RAinc.zeros();
    arma::imat Ainc (nages, incsize); Ainc.zeros();
    arma::imat Einc (nages, incsize); Einc.zeros();
    
    arma::ivec obsInc (data.n_cols); obsInc.zeros();
    arma::ivec Dtempinc (data.n_cols); Dtempinc.zeros();
    arma::vec tempdens (data.n_cols); tempdens.zeros();
    
    // check which output required
    List out (npart * (ndays + 1));
    arma::icube u1_night_t(2, nages, nlads);
    if(saveAll != 0) {
        if(saveAll == 1) {
            for(i = 0; i < npart; i++) {
                u1_night_t.row(0) = u1_night[i].row(6);
                u1_night_t.row(1) = u1_night[i].row(11);
                out[i] = u1_night_t;
            }
        } else {
            for(i = 0; i < npart; i++) {
                out[i] = u1_night[i];
            }
        }
    }
    
    // initialise timer
    Timer timer;
    int timer_cnt = 0;
    double prev_time = 0.0;
    
    // loop over time
    double ll = 0.0;
    for(t = 0; t < ndays; t++) {
        
        // loop over particles
        for(i = 0; i < npart; i++) {
            
            // check for interrupt
            R_CheckUserInterrupt();
            
            // run model and return u1
            discreteStochModel(pars, t - 1, t, &u1_moves, &u1_new[i], &u1_day_new[i], &u1_night_new[i], &N_day, &N_night, &pinf, &origE,  C);
            
            // cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
            //          0,   1,   2,   3,    4,   5,    6,    7,    8,    9,   10,   11
            
            // set weights
            weights(i) = 0.0;
            
            // adjust states according to model discrepancy
            if(MD == 1) {
                // DH (MD on incidence)
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        DHinc(j, l) = u1_new[i](11, j, l) - u1[i](11, j, l);
                        DHinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * DHinc(j, l),
                            a_dis + b_dis * DHinc(j, l),
                            -DHinc(j, l),
                            u1[i](9, j, l) - DHinc(j, l)
                        );
                        if(DHinc(j, l) < 0) Rprintf("DHinc = %d j = %d l = %d\n", DHinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        u1_new[i](11, j, l) = u1[i](11, j, l) + DHinc(j, l);
                        u1_day_new[i](11, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](11, j, (arma::uword) u1_moves(l, 1) - 1) + DHinc(j, l);
                        u1_night_new[i](11, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](11, j, (arma::uword) u1_moves(l, 0) - 1) + DHinc(j, l);
                        if(u1_new[i](11, j, l) < 0) Rprintf("DH = %d j = %d l = %d\n", u1_new[i](11, j, l), j, l);
                        // cu[i](11, j, l) += DHinc(j, l);
                    }
                }
                
                // RH given DH (MD on incidence)
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        RHinc(j, l) = u1_new[i](10, j, l) - u1[i](10, j, l);
                        RHinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * RHinc(j, l),
                            a_dis + b_dis * RHinc(j, l),
                            -RHinc(j, l),
                            u1[i](9, j, l) - DHinc(j, l) - RHinc(j, l)
                        );
                        if(RHinc(j, l) < 0) Rprintf("RHinc = %d j = %d l = %d\n", RHinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        u1_new[i](10, j, l) = u1[i](10, j, l) + RHinc(j, l);
                        u1_day_new[i](10, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](10, j, (arma::uword) u1_moves(l, 1) - 1) + RHinc(j, l);
                        u1_night_new[i](10, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](10, j, (arma::uword) u1_moves(l, 0) - 1) + RHinc(j, l);
                        if(u1_new[i](10, j, l) < 0) Rprintf("RH = %d j = %d l = %d\n", u1_new[i](10, j, l), j, l);
                        // cu[i](10, j, l) += RHinc(j, l);
                    }
                }
                    
                // H given later
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        tempLB = -u1_new[i](9, j, l) + u1[i](9, j, l) - DHinc(j, l) - RHinc(j, l);
                        tempLB = (tempLB > -u1_new[i](9, j, l) ? tempLB:(-u1_new[i](9, j, l)));
                        u1_new[i](9, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](9, j, l),
                            a_dis + b_dis * u1_new[i](9, j, l),
                            tempLB,
                            u1[i](5, j, l) - u1_new[i](9, j, l) + u1[i](9, j, l) - DHinc(j, l) - RHinc(j, l)
                        );
                        Hinc(j, l) = u1_new[i](9, j, l) - u1[i](9, j, l) + DHinc(j, l) + RHinc(j, l);
                        if(Hinc(j, l) < 0) Rprintf("Hinc = %d j = %d l = %d\n", Hinc(j, l), j, l);
                        if(u1_new[i](9, j, l) < 0) Rprintf("H = %d j = %d l = %d\n", u1_new[i](9, j, l), j, l);
                        u1_day_new[i](9, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](9, j, (arma::uword) u1_moves(l, 1) - 1) + Hinc(j, l) - RHinc(j, l) - DHinc(j, l);
                        u1_night_new[i](9, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](9, j, (arma::uword) u1_moves(l, 0) - 1) + Hinc(j, l) - RHinc(j, l) - DHinc(j, l);
                    }
                }
                // // update counts
                // for(j = 0; j < nages; j++) {
                //     for(l = 0; l < incsize; l++) {
                //         cu[i](9, j, l) += Hinc(j, l);
                //     }
                // }
                    
                // DI given H (MD on incidence)
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        DIinc(j, l) = u1_new[i](6, j, l) - u1[i](6, j, l);
                        DIinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * DIinc(j, l),
                            a_dis + b_dis * DIinc(j, l),
                            -DIinc(j, l),
                            u1[i](5, j, l) - Hinc(j, l) - DIinc(j, l)
                        );
                        if(DIinc(j, l) < 0) Rprintf("DIinc = %d j = %d l = %d\n", DIinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        u1_new[i](6, j, l) = u1[i](6, j, l) + DIinc(j, l);
                        u1_day_new[i](6, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](6, j, (arma::uword) u1_moves(l, 1) - 1) + DIinc(j, l);
                        u1_night_new[i](6, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](6, j, (arma::uword) u1_moves(l, 0) - 1) + DIinc(j, l);
                        if(u1_new[i](6, j, l) < 0) Rprintf("DI = %d j = %d l = %d\n", u1_new[i](6, j, l), j, l);
                        // cu[i](6, j, l) += DIinc(j, l);
                    }
                }
                    
                // RI (MD on incidence)
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        RIinc(j, l) = u1_new[i](8, j, l) - u1[i](8, j, l);
                        RIinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * RIinc(j, l),
                            a_dis + b_dis * RIinc(j, l),
                            -RIinc(j, l),
                            u1[i](7, j, l) - RIinc(j, l)
                        );
                        if(RIinc(j, l) < 0) Rprintf("RIinc = %d j = %d l = %d\n", RIinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        u1_new[i](8, j, l) = u1[i](8, j, l) + RIinc(j, l);
                        u1_day_new[i](8, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](8, j, (arma::uword) u1_moves(l, 1) - 1) + RIinc(j, l);
                        u1_night_new[i](8, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](8, j, (arma::uword) u1_moves(l, 0) - 1) + RIinc(j, l);
                        if(u1_new[i](8, j, l) < 0) Rprintf("RI = %d j = %d l = %d\n", u1_new[i](8, j, l), j, l);
                        // cu[i](8, j, l) += RIinc(j, l);
                    }
                }
                    
                // I2 given later
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        tempLB = -u1_new[i](7, j, l) + u1[i](7, j, l) - RIinc(j, l);
                        tempLB = (tempLB > -u1_new[i](7, j, l) ? tempLB:(-u1_new[i](7, j, l)));
                        u1_new[i](7, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](7, j, l),
                            a_dis + b_dis * u1_new[i](7, j, l),
                            tempLB,
                            u1[i](5, j, l) - Hinc(j, l) - DIinc(j, l) - u1_new[i](7, j, l) + u1[i](7, j, l) - RIinc(j, l)
                        );
                        I2inc(j, l) = u1_new[i](7, j, l) - u1[i](7, j, l) + RIinc(j, l);
                        if(I2inc(j, l) < 0) Rprintf("I2inc = %d j = %d l = %d\n", I2inc(j, l), j, l);
                        if(u1_new[i](7, j, l) < 0) Rprintf("I2 = %d j = %d l = %d\n", u1_new[i](7, j, l), j, l);
                        u1_day_new[i](7, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](7, j, (arma::uword) u1_moves(l, 1) - 1) + I2inc(j, l) - RIinc(j, l);
                        u1_night_new[i](7, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](7, j, (arma::uword) u1_moves(l, 0) - 1) + I2inc(j, l) - RIinc(j, l);
                    }
                }
                // update counts
                // for(j = 0; j < nages; j++) {
                //     for(l = 0; l < incsize; l++) {
                //         cu[i](7, j, l) += I2inc(j, l);
                //     }
                // } 
                
                // I1 given later
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        tempLB = -u1_new[i](5, j, l) + u1[i](5, j, l) - I2inc(j, l) - Hinc(j, l) - DIinc(j, l);
                        tempLB = (tempLB > -u1_new[i](5, j, l) ? tempLB:(-u1_new[i](5, j, l)));
                        u1_new[i](5, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](5, j, l),
                            a_dis + b_dis * u1_new[i](5, j, l),
                            tempLB,
                            u1[i](4, j, l) - u1_new[i](5, j, l) + u1[i](5, j, l) - I2inc(j, l) - Hinc(j, l) - DIinc(j, l)
                        );
                        I1inc(j, l) = u1_new[i](5, j, l) - u1[i](5, j, l) + I2inc(j, l) + Hinc(j, l) + DIinc(j, l);
                        if(I1inc(j, l) < 0) Rprintf("I1inc = %d j = %d l = %d\n", I1inc(j, l), j, l);
                        if(u1_new[i](5, j, l) < 0) Rprintf("I1 = %d j = %d l = %d\n", u1_new[i](5, j, l), j, l);
                        u1_day_new[i](5, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](5, j, (arma::uword) u1_moves(l, 1) - 1) + I1inc(j, l) - I2inc(j, l) - DIinc(j, l) - Hinc(j, l);
                        u1_night_new[i](5, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](5, j, (arma::uword) u1_moves(l, 0) - 1) + I1inc(j, l) - I2inc(j, l) - DIinc(j, l) - Hinc(j, l);
                    }
                }
                // update counts
                // for(j = 0; j < nages; j++) {
                //     for(l = 0; l < incsize; l++) {
                //         cu[i](5, j, l) += I1inc(j, l);
                //     }
                // }
                
                // P given later
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        tempLB = -u1_new[i](4, j, l) + u1[i](4, j, l) - I1inc(j, l);
                        tempLB = (tempLB > -u1_new[i](4, j, l) ? tempLB:(-u1_new[i](4, j, l)));
                        u1_new[i](4, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](4, j, l),
                            a_dis + b_dis * u1_new[i](4, j, l),
                            tempLB,
                            u1[i](1, j, l) - u1_new[i](4, j, l) + u1[i](4, j, l) - I1inc(j, l) 
                        );
                        Pinc(j, l) = u1_new[i](4, j, l) - u1[i](4, j, l) + I1inc(j, l);
                        if(Pinc(j, l) < 0) Rprintf("Pinc = %d j = %d l = %d\n", Pinc(j, l), j, l);
                        if(u1_new[i](4, j, l) < 0) Rprintf("P = %d j = %d l = %d\n", u1_new[i](4, j, l), j, l);
                        u1_day_new[i](4, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](4, j, (arma::uword) u1_moves(l, 1) - 1) + Pinc(j, l) - I1inc(j, l);
                        u1_night_new[i](4, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](4, j, (arma::uword) u1_moves(l, 0) - 1) + Pinc(j, l) - I1inc(j, l);
                    }
                }
                // update counts
                // for(j = 0; j < nages; j++) {
                //     for(l = 0; l < incsize; l++) {
                //         cu[i](4, j, l) += Pinc(j, l);
                //     }
                // } 
                
                // RA (MD on incidence)
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        RAinc(j, l) = u1_new[i](3, j, l) - u1[i](3, j, l);
                        RAinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * RAinc(j, l),
                            a_dis + b_dis * RAinc(j, l),
                            -RAinc(j, l),
                            u1[i](2, j, l) - RAinc(j, l)
                        );
                        if(RAinc(j, l) < 0) Rprintf("RAinc = %d j = %d l = %d\n", RAinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        u1_new[i](3, j, l) = u1[i](3, j, l) + RAinc(j, l);
                        u1_day_new[i](3, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](3, j, (arma::uword) u1_moves(l, 1) - 1) + RAinc(j, l);
                        u1_night_new[i](3, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](3, j, (arma::uword) u1_moves(l, 0) - 1) + RAinc(j, l);
                        if(u1_new[i](3, j, l) < 0) Rprintf("RH = %d j = %d l = %d\n", u1_new[i](3, j, l), j, l);
                        // cu[i](3, j, l) += RAinc(j, l);
                    }
                }
                
                // A given later
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        tempLB = -u1_new[i](2, j, l) + u1[i](2, j, l) - RAinc(j, l);
                        tempLB = (tempLB > -u1_new[i](2, j, l) ? tempLB:(-u1_new[i](2, j, l)));
                        u1_new[i](2, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](2, j, l),
                            a_dis + b_dis * u1_new[i](2, j, l),
                            tempLB,
                            u1[i](1, j, l) - Pinc(j, l) - u1_new[i](2, j, l) + u1[i](2, j, l) - RAinc(j, l)
                        );
                        Ainc(j, l) = u1_new[i](2, j, l) - u1[i](2, j, l) + RAinc(j, l);
                        if(Ainc(j, l) < 0) Rprintf("Ainc = %d j = %d l = %d\n", Ainc(j, l), j, l);
                        if(u1_new[i](2, j, l) < 0) Rprintf("A = %d j = %d l = %d\n", u1_new[i](2, j, l), j, l);
                        u1_day_new[i](2, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](2, j, (arma::uword) u1_moves(l, 1) - 1) + Ainc(j, l) - RAinc(j, l);
                        u1_night_new[i](2, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](2, j, (arma::uword) u1_moves(l, 0) - 1) + Ainc(j, l) - RAinc(j, l);
                    }
                }
                // update counts
                // for(j = 0; j < nages; j++) {
                //     for(l = 0; l < incsize; l++) {
                //         cu[i](2, j, l) += Ainc(j, l);
                //     }
                // } 
                
                // E given later
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        tempLB = -u1_new[i](1, j, l) + u1[i](1, j, l) - Ainc(j, l) - Pinc(j, l);
                        tempLB = (tempLB > -u1_new[i](1, j, l) ? tempLB:(-u1_new[i](1, j, l)));
                        u1_new[i](1, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](1, j, l),
                            a_dis + b_dis * u1_new[i](1, j, l),
                            tempLB,
                            u1[i](0, j, l) - u1_new[i](1, j, l) + u1[i](1, j, l) - Ainc(j, l) - Pinc(j, l)
                        );
                        Einc(j, l) = u1_new[i](1, j, l) - u1[i](1, j, l) + Ainc(j, l) + Pinc(j, l);
                        if(Einc(j, l) < 0) Rprintf("Einc = %d j = %d l = %d\n", Einc(j, l), j, l);
                        if(u1_new[i](1, j, l) < 0) Rprintf("E = %d j = %d l = %d\n", u1_new[i](1, j, l), j, l);
                        u1_day_new[i](1, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](1, j, (arma::uword) u1_moves(l, 1) - 1) + Einc(j, l) - Pinc(j, l) - Ainc(j, l);
                        u1_night_new[i](1, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](1, j, (arma::uword) u1_moves(l, 0) - 1) + Einc(j, l) - Pinc(j, l) - Ainc(j, l);
                    }
                }
                // update counts
                // for(j = 0; j < nages; j++) {
                //     for(l = 0; l < incsize; l++) {
                //         cu[i](1, j, l) += Einc(j, l);
                //     }
                // }
                
                // S given later
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < incsize; l++) {
                        u1_new[i](0, j, l) = u1[i](0, j, l) - Einc(j, l);
                        u1_day_new[i](0, j, (arma::uword) u1_moves(l, 1) - 1) = u1_day[i](0, j, (arma::uword) u1_moves(l, 1) - 1) - Einc(j, l);
                        u1_night_new[i](0, j, (arma::uword) u1_moves(l, 0) - 1) = u1_night[i](0, j, (arma::uword) u1_moves(l, 0) - 1) - Einc(j, l);
                        if(u1_new[i](0, j, l) < 0) Rprintf("S = %d j = %d l = %d\n", u1_new[i](0, j, l), j, l);
                        // cu[i](0, j, l) -= Einc(j, l);
                    }
                }
                
                // generate data in correct format for observation error weights
                Dtempinc.zeros();
                for(l = 0; l < incsize; l++) {
                    for(j = 0; j < nages; j++) {
                        Dtempinc((arma::uword) j * nlads + u1_moves(l, 0) - 1) += DIinc(j, l);
                        Dtempinc((arma::uword) nlads * nages + j * nlads + u1_moves(l, 0) - 1) += DHinc(j, l);
                    }
                }
            } else {
                // cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
                //          0,   1,   2,   3,    4,   5,    6,    7,    8,    9,   10,   11
                
                // DH
                DHinc = u1_night_new[i].row(11) - u1_night[i].row(11);
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](11, j, l) += DHinc(j, l);
                // 
                // // RH 
                // RHinc = u1_new[i].row(10) - u1[i].row(10);
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](10, j, l) += RHinc(j, l);
                // 
                // // H
                // Hinc = u1_new[i].row(9) - u1[i].row(9);
                // Hinc = Hinc + DHinc + RHinc;
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](9, j, l) += Hinc(j, l);
                // 
                // DI
                DIinc = u1_night_new[i].row(6) - u1_night[i].row(6);
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](6, j, l) += DIinc(j, l);
                // 
                // // RI
                // RIinc = u1_new[i].row(8) - u1[i].row(8);
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](8, j, l) += RIinc(j, l);
                // 
                // // I2 
                // I2inc = u1_new[i].row(7) - u1[i].row(7);
                // I2inc = I2inc + RIinc;
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](7, j, l) += I2inc(j, l);
                // 
                // // I1 
                // I1inc = u1_new[i].row(5) - u1[i].row(5);
                // I1inc = I1inc + DIinc + Hinc + I2inc;
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](5, j, l) += I1inc(j, l);
                // 
                // // P given later
                // Pinc = u1_new[i].row(4) - u1[i].row(4);
                // Pinc = Pinc + I1inc;
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](4, j, l) += Pinc(j, l);
                // 
                // // RA
                // RAinc = u1_new[i].row(3) - u1[i].row(3);
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](3, j, l) += Ainc(j, l);
                // 
                // // A 
                // Ainc = u1_new[i].row(2) - u1[i].row(2);
                // Ainc = Ainc + RAinc;
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](2, j, l) += RAinc(j, l);
                // 
                // // E
                // Einc = u1_new[i].row(1) - u1[i].row(1);
                // Einc = Einc + Ainc + Pinc;
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](1, j, l) += Einc(j, l);
                // 
                // // S
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) u1_new[i](0, j, l) = u1[i](0, j, l) - Einc(j, l);
                // for(j = 0; j < nages; j++) for(l = 0; l < u1_moves.n_rows; l++) cu[i](0, j, l) -= Einc(j, l);
                
                // generate data in correct format for observation error weights
                Dtempinc.zeros();
                for(l = 0; l < nlads; l++) {
                    for(j = 0; j < nages; j++) {
                        Dtempinc((arma::uword) j * nlads + l) += DIinc(j, l);
                        Dtempinc((arma::uword) nlads * nages + j * nlads + l) += DHinc(j, l);
                    }
                }
            }
            
            // calculate log observation error weights
            obsInc = data.row(t).t();
            for(l = 0; l < Dtempinc.n_elem; l++) {
                tempdens(l) = ldtskellam_cpp(
                    obsInc(l) - Dtempinc(l),
                    a1 + b * Dtempinc(l),
                    a2 + b * Dtempinc(l),
                    -Dtempinc(l),
                    obsInc(l)
                );
            }
            weights(i) += log_sum_exp(tempdens, 1);
        }
        
        // calculate log-likelihood contribution
        ll += log_sum_exp(weights, 1);
        
        // // if zero likelihood then return
        // if(!is.finite(ll)) {
        //     if(saveAll == 0) {
        //         return List::create(Named("ll") = ll);
        //     } else {
        //         return List::create(Named("ll") = ll, _["particles"] = out);
        //     }
        // }
        
        // normalise weights
        wnorm = log_sum_exp(weights, 0);
        weights = exp(weights - wnorm);
        
        // resample
        rmultinom(npart, weights.begin(), npart, inds.begin());
        l = 0;
        for(i = 0; i < npart; i++) {
            j = 0;
            while(j < inds(i)) {
                u1[l] = u1_new[i];
                u1_day[l] = u1_day_new[i];
                u1_night[l] = u1_night_new[i];
                l++;
                j++;
            }
        }
        if(l != npart) Rprintf("Something wrong with re-sampling step\n");
        
        // copy in order to pass by reference
        for(i = 0; i < npart; i++) {
            u1_new[i] = u1[i];
            u1_day_new[i] = u1_day[i];
            u1_night_new[i] = u1_night[i];
        }
        
        // save particles if necessary
        if(saveAll != 0) {
            if(saveAll == 1) {
                for(i = 0; i < npart; i++) {
                    u1_night_t.row(0) = u1_night[i].row(6);
                    u1_night_t.row(1) = u1_night[i].row(11);
                    out[i + npart * (t + 1)] = u1_night_t;
                }
            } else {
                for(i = 0; i < npart; i++) {
                    out[i + npart * (t + 1)] = u1_night[i];
                }
            }
        }
        
        //calculate block run time
        timer.step("");
        NumericVector res(timer);
        
        Rprintf("t = %d / %d time = %.2f secs \n", t, ndays, (res[timer_cnt] / 1e9) - prev_time);
        
        //reset timer and acceptance rate counter
        prev_time = res[timer_cnt] / 1e9;
        timer_cnt++;
    }
    if(saveAll == 0) {
        return List::create(Named("ll") = ll);
    } else {
        return List::create(Named("ll") = ll, _["particles"] = out);
    }
}
