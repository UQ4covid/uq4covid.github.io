// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::plugins(openmp)]]

// [[Rcpp::depends(sitmo)]]

#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include <sitmo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// log-sum-exp function to prevent numerical overflow
double log_sum_exp(arma::vec x, int mn = 0) {
    double maxx = max(x);
    double y = maxx + log(sum(exp(x - maxx)));
    if(mn == 1) y = y - log(x.n_elem);
    return y;
}

// Poisson RNG using inverse transform method
// (to try to circumvent non thread-safe RNG in R)
int rpois_cpp (double lambda, sitmo::prng &eng) {
    if(lambda < 0.0) {
        stop("'lambda' must be > 0 in rpois\n");
    }
    double mx = sitmo::prng::max();
    double u = log(eng()) - log(mx);
    double temp = -lambda;
    double temp1 = 0.0;
    double maxt = 0.0;
    double lk = 0.0;
    int k = 0;
    while(temp < u) {
        k++;
        temp1 = k * log(lambda) - lambda;
        lk += log(k);
        temp1 -= lk;
        maxt = (temp > temp1 ? temp:temp1);
        temp = maxt + log(exp(temp - maxt) + exp(temp1 - maxt));
    }
    return k;
}

// Binomial RNG using inverse transform method
// (to try to circumvent non thread-safe RNG in R)
int rbinom_cpp (int n, double p, sitmo::prng &eng) {
    if(n < 0) {
        Rprintf("n = %d\n", n);
        stop("'n' must be >= 0 in rbinom\n");
    }
    if(p < 0.0 || p > 1.0) {
        Rprintf("p = %f\n", p);
        stop("Must have 0 <= p <= 1 in rbinom\n");
    }
    if(n == 0) return 0;
    double mx = sitmo::prng::max();
    double u = log(eng()) - log(mx);
    int k;
    double ln = 0.0, lx = 0.0, lnmx = 0.0;
    for(k = 1; k <= n; k++) {
        ln += log(k);
    }
    lnmx = ln;
    double temp = ln - lx - lnmx + n * log(1.0 - p);
    double temp1 = 0.0, maxt = 0.0;
    k = 0;
    while(temp < u) {
        k++;
        lx += log(k);
        lnmx -= log(n - (k - 1));
        temp1 = ln - lx - lnmx + k * log(p) + (n - k) * log(1.0 - p);
        maxt = (temp > temp1 ? temp:temp1);
        temp = maxt + log(exp(temp - maxt) + exp(temp1 - maxt));
    }
    return k;
}

// Multinomial RNG for n = 1 using inverse transform method
// (to try to circumvent non thread-safe RNG in R)
int rmultinom_cpp (arma::vec p, sitmo::prng &eng) {
    if(abs(sum(p) - 1.0) > 1e-15) {
        Rprintf("sum(p) = %e\n", sum(p));
        stop("Must have sum(p) == 1 in rmultinom\n");
    }
    if(any(p) > 1.0 || any(p) < 0.0) {
        stop("Must have all 0 <= p <= 1 in rmultinom\n");
    }
    double mx = sitmo::prng::max();
    double u = eng() / mx;
    double temp = p(0);
    int k = 0;
    while(temp < u) {
        k++;
        temp += p(k);
    }
    return k;
}  

// log Skellam CDF (note, no checks on inputs)
// adapted from "pskellam" source code by Patrick Brown (mistakes are mine)
double lpskellam_cpp(int x, double lambda1, double lambda2) {
    // different formulas for negative & nonnegative x (zero lambda is OK)
    double ldens = 0.0;
    if(x < 0) {
        ldens = R::pnchisq(2.0 * lambda2, -2.0 * x, 2.0 * lambda1, 1, 1);
    } else {
        ldens = R::pnchisq(2.0 * lambda1, 2.0 * (x + 1), 2.0 * lambda2, 0, 1);
    }
    if(!arma::is_finite(ldens)) {
        Rprintf("Non-finite pdensity in lpskellam = %f x = %d lambda1 = %f lambda2 = %f\n", ldens, x, lambda1, lambda2);
    }
    return ldens;
}

// log truncated Skellam density (note, no checks on inputs)
// adapted from "dskellam" source code by Patrick Brown (mistakes are mine)
double ldtskellam_cpp(int x, double lambda1, double lambda2, char *str, int LB = 0, int UB = -1, int print = 0) {
    
    // if UB < LB then returns untruncated density
    
    // check for this instance only
    if(LB > UB) {
        stop("ldt LB > UB\n");
    }
    
    if(print == 1) Rprintf("State: %s\n", str);
    
    // create array for Bessel (trying to deal with recursive
    // gc error message that may be coming from memory allocation
    // in bessel_i - hence swapping to bessel_i_ex and controlling
    // vector allocation directly here through calloc and free -
    // OK to use std::calloc here because object not going
    // back to R I think)
    double *bi = (double *) calloc (floor(abs(x)) + 1, sizeof(double));
    
    // from "skellam" source code: log(besselI(y, nu)) == y + log(besselI(y, nu, TRUE))
    double ldens = -(lambda1 + lambda2) + (x / 2.0) * (log(lambda1) - log(lambda2)) +
        log(R::bessel_i_ex(2.0 * sqrt(lambda1 * lambda2), abs(x), 2, bi)) + 2.0 * sqrt(lambda1 * lambda2);
    
    // free memory from the heap
    free(bi);
    
    // if non-finite density, then try difference of CDFs
    if(!arma::is_finite(ldens)) {
        Rprintf("Trying difference of CDFs in ldtskellam\n");
        arma::vec norm(2); norm.zeros();
        norm(0) = lpskellam_cpp(x - 1, lambda1, lambda2);
        norm(1) = lpskellam_cpp(x, lambda1, lambda2);                
        ldens = norm(1) + log(1.0 - exp(norm(0) - norm(1)));
    }
        
    // if truncated then adjust log-density
    if(UB >= LB) {
        if(x < LB || x > UB) {
            Rprintf("'x' must be in bounds in dtskellam_cpp: x = %d LB = %d UB = %d\n", x, LB, UB);
            stop("");
        }
        if(arma::is_finite(ldens)) {
            // declare variables
            int normsize = UB - LB + 1;
            if(normsize > 1) {
                arma::vec norm(2); norm.zeros();
                norm(0) = lpskellam_cpp(LB - 1, lambda1, lambda2);
                norm(1) = lpskellam_cpp(UB, lambda1, lambda2);
                //Rprintf("norm(0) = %f norm(1) = %f\n", norm(0), norm(1));               
                double lnorm = norm(1) + log(1.0 - exp(norm(0) - norm(1)));
                ldens -= lnorm;
            } else {
                ldens = 0.0;
            }
        }
    }
    if(!arma::is_finite(ldens)) {
        Rprintf("Non-finite density in ldtskellam = %f x = %d l1 = %f l2 = %f LB = %d UB = %d\n", ldens, x, lambda1, lambda2, LB, UB);
    }
    return ldens;
}

// truncated Skellam sampler
int rtskellam_cpp(double lambda1, double lambda2, char *str, sitmo::prng &eng, int LB = 0, int UB = -1) {
    
    // if UB < LB then returns untruncated density
    
    // check in this setting:
    if(LB > UB) {
        stop("rdt LB > UB\n");
    }
    
    // declare variables
    int x = 0;
    double mx = sitmo::prng::max();
    
    // if bounds are the same, then return the
    // only viable value
    if(LB == UB) {
        x = LB;
    } else {
        // draw untruncated sample
        x = rpois_cpp(lambda1, eng) - rpois_cpp(lambda2, eng);
        
        // rejection sample if necessary
        if(UB > LB) {
            int k = 0;
            int ntries = 1000;
            while((x < LB || x > UB) && k < ntries) {
                x = rpois_cpp(lambda1, eng) - rpois_cpp(lambda2, eng);
                k++;
            }
            // if rejection sampling doesn't work
            // then try inverse transform sampling
            if(k == ntries) {
                // Rprintf("LB = %d UB = %d\n", LB, UB);
                double u = eng() / mx;
                k = 0;
                double xdens = exp(ldtskellam_cpp(LB + k, lambda1, lambda2, str, LB, UB, 1));
                while(u > xdens && k < (UB - LB)) {
                    k++;
                    xdens += exp(ldtskellam_cpp(LB + k, lambda1, lambda2, str, LB, UB, 1));
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
void discreteStochModel(int ipart, arma::vec &pars, int tstart, int tstop, 
                        arma::imat &u1_moves, std::vector<arma::icube> &u1, arma::icube &u1_day, arma::icube &u1_night,
                        arma::imat &N_day, arma::imat &N_night, arma::mat &pinf, arma::imat &origE, arma::mat &C, sitmo::prng &eng) {
    
    // u1_moves is a matrix with columns: LADfrom, LADto
    // u1 is a 3D array with dimensions: nclasses x nages x nmoves
    //          each row of u1 must match u1_moves
    // u1_day/night are 3D arrays with dimensions: nclasses x nages x nlads
    // N_day/night are nlad x nage matrices of population counts
    // pinf is nages x nlads auxiliary matrix
    // origE is nages x nmoves auxiliary matrix
    
    // set up auxiliary matrix for counts
    arma::uword nclasses = (arma::uword) u1[ipart].n_rows;
    arma::uword nages = (arma::uword) u1[ipart].n_cols;
    int k, n;
    arma::uword i, j, l;
    
    // reconstruct day/night counts
    u1_day.zeros();
    u1_night.zeros();
    for(i = 0; i < u1_moves.n_rows; i++) {
        for(j = 0; j < nages; j++) {
            for(l = 0; l < nclasses; l++) {
                u1_day(l, j, (arma::uword) u1_moves(i, 1) - 1) += u1[ipart](l, j, i);
                u1_night(l, j, (arma::uword) u1_moves(i, 0) - 1) += u1[ipart](l, j, i);
            }
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
    
    int tcurr = 0;
    
    // set up vector of number of infectives
    arma::vec uinf(nages);
    
    // set up transmission rates
    arma::mat beta(nages, 1);
    
    // set up auxiliary variables
    arma::vec mprobsE(3);
    arma::vec mprobsI1(4);
    arma::vec mprobsH(3);
    tcurr++;
    tstart++;
    
    while(tstart <= tstop) {
        
        // classes are: S, E, A, RA, P, I1, DI, I2, RI, H, RH, DH
        //              0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
        
        // save current number of infectives for later transitions
        for(i = 0; i < u1_moves.n_rows; i++) {
            for(j = 0; j < nages; j++) {
                origE(j, i) = u1[ipart](1, j, i);
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
                pinf(j, i) = 1.0 - exp(-beta(j, 0));
                pinf(j, i) = (pinf(j, i) < 0.0 ? 0.0:pinf(j, i));
                pinf(j, i) = (pinf(j, i) > 1.0 ? 1.0:pinf(j, i));
            }
        }
        // transmission events (day), loop over network
        for(i = 0; i < u1_moves.n_rows; i++) {
            for(j = 0; j < nages; j++) {
                k = rbinom_cpp(u1[ipart](0, j, i), pinf(j, (arma::uword) u1_moves(i, 1) - 1), eng);
                u1[ipart](0, j, i) -= k;
                u1[ipart](1, j, i) += k;
                u1_day(0, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(1, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(0, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(1, j, (arma::uword) u1_moves(i, 0) - 1) += k;
            }
        }
        
        // transmission probabilities (night), loop over LADs
        for(i = 0; i < u1_night.n_slices; i++) {
            
            // update infective counts for rate
            for(j = 0; j < nages; j++) {
                uinf(j) = (double) nuA * u1_night(2, j, i) + nu * (u1_night(4, j, i) + u1_night(5, j, i) + u1_night(7, j, i));
            }
            
            // SE
            beta = 0.3 * C * (uinf / N_night.col(i));
            for(j = 0; j < nages; j++) {
                pinf(j, i) = 1.0 - exp(-beta(j, 0));
                pinf(j, i) = (pinf(j, i) < 0.0 ? 0.0:pinf(j, i));
                pinf(j, i) = (pinf(j, i) > 1.0 ? 1.0:pinf(j, i));
            }
        }
        // transmission events (night), loop over network
        for(i = 0; i < u1_moves.n_rows; i++) {
            for(j = 0; j < nages; j++) {
                k = rbinom_cpp(u1[ipart](0, j, i), pinf(j, (arma::uword) u1_moves(i, 0) - 1), eng);
                u1[ipart](0, j, i) -= k;
                u1[ipart](1, j, i) += k;
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
            for(i = 0; i < u1_moves.n_rows; i++) {
                
                // H out
                n = u1[ipart](9, j, i);
                if(n > 0) {
                    for(l = 0; l < n; l++) {
                        k = rmultinom_cpp(mprobsH, eng);
                        u1[ipart](9, j, i) -= (k < 2 ? 1:0);
                        u1[ipart](10, j, i) += (k == 0 ? 1:0);
                        u1[ipart](11, j, i) += (k == 1 ? 1:0);
                        u1_day(9, j, (arma::uword) u1_moves(i, 1) - 1) -= (k < 2 ? 1:0);
                        u1_day(10, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 0 ? 1:0);
                        u1_day(11, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 1 ? 1:0);
                        u1_night(9, j, (arma::uword) u1_moves(i, 0) - 1) -= (k < 2 ? 1:0);
                        u1_night(10, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 0 ? 1:0);
                        u1_night(11, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 1 ? 1:0);
                    }
                }
                
                // I2RI
                k = rbinom_cpp(u1[ipart](7, j, i), probI2(j), eng);
                u1[ipart](7, j, i) -= k;
                u1[ipart](8, j, i) += k;
                u1_day(7, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(8, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(7, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(8, j, (arma::uword) u1_moves(i, 0) - 1) += k;
                
                // I1 out
                n = u1[ipart](5, j, i);
                if(n > 0) {
                    for(l = 0; l < n; l++) {
                        k = rmultinom_cpp(mprobsI1, eng);
                        u1[ipart](5, j, i) -= (k < 3 ? 1:0);
                        u1[ipart](9, j, i) += (k == 0 ? 1:0);
                        u1[ipart](7, j, i) += (k == 1 ? 1:0);
                        u1[ipart](6, j, i) += (k == 2 ? 1:0);
                        u1_day(5, j, (arma::uword) u1_moves(i, 1) - 1) -= (k < 3 ? 1:0);
                        u1_day(9, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 0 ? 1:0);
                        u1_day(7, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 1 ? 1:0);
                        u1_day(6, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 2 ? 1:0);
                        u1_night(5, j, (arma::uword) u1_moves(i, 0) - 1) -= (k < 3 ? 1:0);
                        u1_night(9, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 0 ? 1:0);
                        u1_night(7, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 1 ? 1:0);
                        u1_night(6, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 2 ? 1:0);
                    }
                }
                
                // PI1
                k = rbinom_cpp(u1[ipart](4, j, i), probP(j), eng);
                u1[ipart](4, j, i) -= k;
                u1[ipart](5, j, i) += k;
                u1_day(4, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(5, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(4, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(5, j, (arma::uword) u1_moves(i, 0) - 1) += k;
                
                // ARA
                k = rbinom_cpp(u1[ipart](2, j, i), probA(j), eng);
                u1[ipart](2, j, i) -= k;
                u1[ipart](3, j, i) += k;
                u1_day(2, j, (arma::uword) u1_moves(i, 1) - 1) -= k;
                u1_day(3, j, (arma::uword) u1_moves(i, 1) - 1) += k;
                u1_night(2, j, (arma::uword) u1_moves(i, 0) - 1) -= k;
                u1_night(3, j, (arma::uword) u1_moves(i, 0) - 1) += k;
                
                // E out
                if(origE(j, i) > 0) {
                    for(l = 0; l < origE(j, i); l++) {
                        k = rmultinom_cpp(mprobsE, eng);
                        u1[ipart](1, j, i) -= (k < 2 ? 1:0);
                        u1[ipart](2, j, i) += (k == 0 ? 1:0);
                        u1[ipart](4, j, i) += (k == 1 ? 1:0);
                        u1_day(1, j, (arma::uword) u1_moves(i, 1) - 1) -= (k < 2 ? 1:0);
                        u1_day(2, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 0 ? 1:0);
                        u1_day(4, j, (arma::uword) u1_moves(i, 1) - 1) += (k == 1 ? 1:0);
                        u1_night(1, j, (arma::uword) u1_moves(i, 0) - 1) -= (k < 2 ? 1:0);
                        u1_night(2, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 0 ? 1:0);
                        u1_night(4, j, (arma::uword) u1_moves(i, 0) - 1) += (k == 1 ? 1:0);
                    }
                }
            }
        }
        
        // update time 
        tcurr++;
        tstart++;
    }
    return;
}

//// [[Rcpp::export]]
//void testsitmo () {
//    // check RNGs
//    int ncores = 2;
//    arma::vec seeds(ncores);
//    for(int i = 0; i < ncores; i++) {
//        seeds(i) = R::rnorm(0.0, 100.0);
//    }
//    double mx = sitmo::prng::max();
//    Rprintf("Without adjusting seed state:\n");
//    for(int j = 0; j < 2; j++) {
//#pragma omp parallel for default(none) shared(j, ncores, seeds, mx)
//        for(int i = 0; i < ncores; i++) {
//            uint32_t coreseed = static_cast<uint32_t>(seeds((arma::uword) i));
//            sitmo::prng eng(coreseed);
//            Rprintf("j = %d i = %d runif = %f\n", j, i, eng() / mx);
//        }
//    }
//    Rprintf("With adjusting seed state:\n");
//    for(int j = 0; j < 2; j++) {
//#pragma omp parallel for default(none) shared(j, ncores, seeds, mx)
//        for(int i = 0; i < ncores; i++) {
//            uint32_t coreseed = static_cast<uint32_t>(seeds((arma::uword) i));
//            sitmo::prng eng(coreseed);
//            Rprintf("j = %d i = %d runif = %f\n", j, i, eng() / mx);
//            seeds((arma::uword) i) = eng();
//        }
//    }
//    return;
//}

//// [[Rcpp::export]]
//List testRNGs () {
//    // check RNGs
//    uint32_t coreseed = static_cast<uint32_t>(R::rnorm(0.0, 100.0));
//    sitmo::prng eng(coreseed);
//    arma::ivec resBin (10000);
//    arma::ivec resPois (10000);
//    arma::ivec resMulti (10000);
//    arma::vec probs = {0.1, 0.2, 0.05, 0.35, 0.05, 0.05, 0.1, 0.1};
//    for(arma::uword i = 0; i < 10000; i++) {
//        resBin(i) = rbinom_cpp(30, 0.89, eng);
//        resPois(i) = rpois_cpp(0.8, eng);
//        resMulti(i) = rmultinom_cpp(probs, eng);
//    }
//    return List::create(Named("bin") = resBin, _["pois"] = resPois, _["multi"] = resMulti);
//}

// [[Rcpp::export]]
List PF_cpp (arma::vec pars, arma::mat C, arma::imat data, arma::uword nclasses, arma::uword nages, arma::uword nlads, arma::imat u1_moves, 
         arma::icube u1_comb, arma::uword ndays, arma::uword npart, int MD, double a1, double a2, double b, double a_dis, 
         double b_dis, int saveAll) {
    
    // set counters
    arma::uword i, j, l, k, t;
    int tempLB = 0;
    
    // split u1 up into different LADs
    std::vector<arma::icube> u1(npart);
    std::vector<arma::icube> u1_new(npart);
    arma::icube u1_night_full(nclasses, nages, nlads); u1_night_full.zeros();
    arma::icube u1_night_reduced(2, nages, nlads); u1_night_reduced.zeros();
    arma::imat N_day(nages, nlads); N_day.zeros();
    arma::imat N_night(nages, nlads); N_night.zeros();
    for(i = 0; i < u1_moves.n_rows; i++) {
        for(j = 0; j < nages; j++) {
            for(l = 0; l < nclasses; l++) {
                N_day(j, (arma::uword) u1_moves(i, 1) - 1) += u1_comb(l, j, i);
                N_night(j, (arma::uword) u1_moves(i, 0) - 1) += u1_comb(l, j, i);
            }
        }
    }
    for(i = 0; i < npart; i++) {
        u1[i] = u1_comb;
        u1_new[i] = u1_comb;
    }
    
    // set up weight vector
    arma::vec weights (npart);
    arma::ivec inds(npart);
    double wnorm = 0.0;
    
    // set up auxiliary objects    
    arma::ivec obsInc (data.n_cols); obsInc.zeros();
    
    // check which output required
    List out (npart * (ndays + 1));
    if(saveAll != 0) {
        if(saveAll == 1) {
            for(i = 0; i < npart; i++) {
                // extract just counts for DI and DH
                u1_night_reduced.zeros();
                for(l = 0; l < u1_moves.n_rows; l++) {
                    for(j = 0; j < nages; j++) {
                        u1_night_reduced(0, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](6, j, l);
                        u1_night_reduced(1, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](11, j, l);
                    }
                }
                out[i] = u1_night_reduced;
            }
        } else {
            for(i = 0; i < npart; i++) {
                u1_night_full.zeros();
                for(l = 0; l < u1_moves.n_rows; l++) {
                    for(j = 0; j < nages; j++) {
                        for(k = 0; k < nclasses; k++) {
                            u1_night_full(k, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](k, j, l);
                        }
                    }
                }
                out[i] = u1_night_full;
            }
        }
    }
    
    // initialise timer
    Timer timer;
    int timer_cnt = 0;
    double prev_time = 0.0;
    
    // sample seeds to set up thread-safe PRNGs
    int ncores = 1;
#ifdef _OPENMP
    ncores = omp_get_max_threads();
    omp_set_num_threads(ncores);
#endif
    arma::vec seeds(ncores);
    for(i = 0; i < ncores; i++) {
        seeds(i) = R::rnorm(0.0, 100.0);
    }
    uint32_t coreseedSerial = static_cast<uint32_t>(R::rnorm(0.0, 100.0));
    sitmo::prng engSerial(coreseedSerial);
    
    // loop over time
    double ll = 0.0;
    for(t = 0; t < ndays; t++) {
            
        // check for interrupt
        R_CheckUserInterrupt();
        
        // extract data
        obsInc = data.row(t).t();
        
        // loop over particles
#ifdef _OPENMP
#pragma omp parallel for default(none) private(j, l, tempLB) shared(seeds, npart, u1_moves, nages, nclasses, nlads, data, C, N_night, N_day, u1, u1_new, t, pars, weights, MD, a_dis, b_dis, a1, a2, b, obsInc)
#endif
        for(i = 0; i < npart; i++) {
    
            // set up print string for debugging
            char str1[80];
	        std::strcpy(str1, "setup");
 
            // set up thread-safe RNG
            uint32_t coreseed = static_cast<uint32_t>(seeds(0));
#ifdef _OPENMP
            coreseed = static_cast<uint32_t>(seeds((arma::uword) omp_get_thread_num()));
#endif
            sitmo::prng eng(coreseed);
        
            // set up auxiliary objects
            arma::imat DHinc (nages, u1_moves.n_rows); DHinc.zeros();
            arma::imat DIinc (nages, u1_moves.n_rows); DIinc.zeros();
            
            arma::ivec Dtempinc (data.n_cols); Dtempinc.zeros();
            arma::vec tempdens (data.n_cols); tempdens.zeros();
            
            arma::mat pinf(nages, nlads); pinf.zeros();
            arma::imat origE(nages, u1_moves.n_rows); origE.zeros();
            arma::icube u1_day(nclasses, nages, nlads); u1_day.zeros();
            arma::icube u1_night(nclasses, nages, nlads); u1_night.zeros();
            
            // run model and return u1
            discreteStochModel((int) i, pars, t - 1, t, u1_moves, u1_new, u1_day, u1_night, N_day, N_night, pinf, origE, C, eng);
            
            // cols: c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
            //          0,   1,   2,   3,    4,   5,    6,    7,    8,    9,   10,   11
            
            // set weights
            weights(i) = 0.0;
            
            // adjust states according to model discrepancy
            if(MD == 1) {
                // create auxiliary objects
                arma::imat RHinc (nages, u1_moves.n_rows); RHinc.zeros();
                arma::imat Hinc (nages, u1_moves.n_rows); Hinc.zeros();
                arma::imat RIinc (nages, u1_moves.n_rows); RIinc.zeros();
                arma::imat I2inc (nages, u1_moves.n_rows); I2inc.zeros();
                arma::imat I1inc (nages, u1_moves.n_rows); I1inc.zeros();
                arma::imat Pinc (nages, u1_moves.n_rows); Pinc.zeros();
                arma::imat RAinc (nages, u1_moves.n_rows); RAinc.zeros();
                arma::imat Ainc (nages, u1_moves.n_rows); Ainc.zeros();
                arma::imat Einc (nages, u1_moves.n_rows); Einc.zeros();
                
                // DH (MD on incidence)
                std::strcpy(str1, "DHinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        DHinc(j, l) = u1_new[i](11, j, l) - u1[i](11, j, l);
                        DHinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * DHinc(j, l),
                            a_dis + b_dis * DHinc(j, l),
                            str1,
                            eng,
                            -DHinc(j, l),
                            u1[i](9, j, l) - DHinc(j, l)
                        );
                        if(DHinc(j, l) < 0) Rprintf("DHinc = %d j = %d l = %d\n", DHinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](11, j, l) = u1[i](11, j, l) + DHinc(j, l);
                        if(u1_new[i](11, j, l) < 0) Rprintf("DH = %d j = %d l = %d\n", u1_new[i](11, j, l), j, l);
                    }
                }
                
                // RH given DH (MD on incidence)
                std::strcpy(str1, "RHinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        RHinc(j, l) = u1_new[i](10, j, l) - u1[i](10, j, l);
                        RHinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * RHinc(j, l),
                            a_dis + b_dis * RHinc(j, l),
                            str1,
                            eng,
                            -RHinc(j, l),
                            u1[i](9, j, l) - DHinc(j, l) - RHinc(j, l)
                        );
                        if(RHinc(j, l) < 0) Rprintf("RHinc = %d j = %d l = %d\n", RHinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](10, j, l) = u1[i](10, j, l) + RHinc(j, l);
                        if(u1_new[i](10, j, l) < 0) Rprintf("RH = %d j = %d l = %d\n", u1_new[i](10, j, l), j, l);
                    }
                }
                    
                // H given later
                std::strcpy(str1, "Hinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](9, j, l) + u1[i](9, j, l) - DHinc(j, l) - RHinc(j, l);
                        tempLB = (tempLB > -u1_new[i](9, j, l) ? tempLB:(-u1_new[i](9, j, l)));
                        u1_new[i](9, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](9, j, l),
                            a_dis + b_dis * u1_new[i](9, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](5, j, l) - u1_new[i](9, j, l) + u1[i](9, j, l) - DHinc(j, l) - RHinc(j, l)
                        );
                        Hinc(j, l) = u1_new[i](9, j, l) - u1[i](9, j, l) + DHinc(j, l) + RHinc(j, l);
                        if(Hinc(j, l) < 0) Rprintf("Hinc = %d j = %d l = %d\n", Hinc(j, l), j, l);
                        if(u1_new[i](9, j, l) < 0) Rprintf("H = %d j = %d l = %d\n", u1_new[i](9, j, l), j, l);
                    }
                }
                    
                // DI given H (MD on incidence)
                std::strcpy(str1, "DIinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        DIinc(j, l) = u1_new[i](6, j, l) - u1[i](6, j, l);
                        DIinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * DIinc(j, l),
                            a_dis + b_dis * DIinc(j, l),
                            str1,
                            eng,
                            -DIinc(j, l),
                            u1[i](5, j, l) - Hinc(j, l) - DIinc(j, l)
                        );
                        if(DIinc(j, l) < 0) Rprintf("DIinc = %d j = %d l = %d\n", DIinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](6, j, l) = u1[i](6, j, l) + DIinc(j, l);
                        if(u1_new[i](6, j, l) < 0) Rprintf("DI = %d j = %d l = %d\n", u1_new[i](6, j, l), j, l);
                    }
                }
                    
                // RI (MD on incidence)
                std::strcpy(str1, "RIinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        RIinc(j, l) = u1_new[i](8, j, l) - u1[i](8, j, l);
                        RIinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * RIinc(j, l),
                            a_dis + b_dis * RIinc(j, l),
                            str1,
                            eng,
                            -RIinc(j, l),
                            u1[i](7, j, l) - RIinc(j, l)
                        );
                        if(RIinc(j, l) < 0) Rprintf("RIinc = %d j = %d l = %d\n", RIinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](8, j, l) = u1[i](8, j, l) + RIinc(j, l);
                        if(u1_new[i](8, j, l) < 0) Rprintf("RI = %d j = %d l = %d\n", u1_new[i](8, j, l), j, l);
                    }
                }
                    
                // I2 given later
                std::strcpy(str1, "I2inc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](7, j, l) + u1[i](7, j, l) - RIinc(j, l);
                        tempLB = (tempLB > -u1_new[i](7, j, l) ? tempLB:(-u1_new[i](7, j, l)));
                        u1_new[i](7, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](7, j, l),
                            a_dis + b_dis * u1_new[i](7, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](5, j, l) - Hinc(j, l) - DIinc(j, l) - u1_new[i](7, j, l) + u1[i](7, j, l) - RIinc(j, l)
                        );
                        I2inc(j, l) = u1_new[i](7, j, l) - u1[i](7, j, l) + RIinc(j, l);
                        if(I2inc(j, l) < 0) Rprintf("I2inc = %d j = %d l = %d\n", I2inc(j, l), j, l);
                        if(u1_new[i](7, j, l) < 0) Rprintf("I2 = %d j = %d l = %d\n", u1_new[i](7, j, l), j, l);
                    }
                }
                
                // I1 given later
                std::strcpy(str1, "I1inc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](5, j, l) + u1[i](5, j, l) - I2inc(j, l) - Hinc(j, l) - DIinc(j, l);
                        tempLB = (tempLB > -u1_new[i](5, j, l) ? tempLB:(-u1_new[i](5, j, l)));
                        u1_new[i](5, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](5, j, l),
                            a_dis + b_dis * u1_new[i](5, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](4, j, l) - u1_new[i](5, j, l) + u1[i](5, j, l) - I2inc(j, l) - Hinc(j, l) - DIinc(j, l)
                        );
                        I1inc(j, l) = u1_new[i](5, j, l) - u1[i](5, j, l) + I2inc(j, l) + Hinc(j, l) + DIinc(j, l);
                        if(I1inc(j, l) < 0) Rprintf("I1inc = %d j = %d l = %d\n", I1inc(j, l), j, l);
                        if(u1_new[i](5, j, l) < 0) Rprintf("I1 = %d j = %d l = %d\n", u1_new[i](5, j, l), j, l);
                    }
                }
                
                // P given later
                std::strcpy(str1, "Pinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](4, j, l) + u1[i](4, j, l) - I1inc(j, l);
                        tempLB = (tempLB > -u1_new[i](4, j, l) ? tempLB:(-u1_new[i](4, j, l)));
                        u1_new[i](4, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](4, j, l),
                            a_dis + b_dis * u1_new[i](4, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](1, j, l) - u1_new[i](4, j, l) + u1[i](4, j, l) - I1inc(j, l)
                        );
                        Pinc(j, l) = u1_new[i](4, j, l) - u1[i](4, j, l) + I1inc(j, l);
                        if(Pinc(j, l) < 0) Rprintf("Pinc = %d j = %d l = %d\n", Pinc(j, l), j, l);
                        if(u1_new[i](4, j, l) < 0) Rprintf("P = %d j = %d l = %d\n", u1_new[i](4, j, l), j, l);
                    }
                }
                
                // RA (MD on incidence)
                std::strcpy(str1, "RAinc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        RAinc(j, l) = u1_new[i](3, j, l) - u1[i](3, j, l);
                        RAinc(j, l) += rtskellam_cpp(
                            a_dis + b_dis * RAinc(j, l),
                            a_dis + b_dis * RAinc(j, l),
                            str1,
                            eng,
                            -RAinc(j, l),
                            u1[i](2, j, l) - RAinc(j, l)
                        );
                        if(RAinc(j, l) < 0) Rprintf("RAinc = %d j = %d l = %d\n", RAinc(j, l), j, l);
                    }
                }
                // update counts
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](3, j, l) = u1[i](3, j, l) + RAinc(j, l);
                        if(u1_new[i](3, j, l) < 0) Rprintf("RH = %d j = %d l = %d\n", u1_new[i](3, j, l), j, l);
                    }
                }
                
                // A given later
                std::strcpy(str1, "Ainc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](2, j, l) + u1[i](2, j, l) - RAinc(j, l);
                        tempLB = (tempLB > -u1_new[i](2, j, l) ? tempLB:(-u1_new[i](2, j, l)));
                        u1_new[i](2, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](2, j, l),
                            a_dis + b_dis * u1_new[i](2, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](1, j, l) - Pinc(j, l) - u1_new[i](2, j, l) + u1[i](2, j, l) - RAinc(j, l)
                        );
                        Ainc(j, l) = u1_new[i](2, j, l) - u1[i](2, j, l) + RAinc(j, l);
                        if(Ainc(j, l) < 0) Rprintf("Ainc = %d j = %d l = %d\n", Ainc(j, l), j, l);
                        if(u1_new[i](2, j, l) < 0) Rprintf("A = %d j = %d l = %d\n", u1_new[i](2, j, l), j, l);
                    }
                }
                
                // E given later
                std::strcpy(str1, "Einc");
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        tempLB = -u1_new[i](1, j, l) + u1[i](1, j, l) - Ainc(j, l) - Pinc(j, l);
                        tempLB = (tempLB > -u1_new[i](1, j, l) ? tempLB:(-u1_new[i](1, j, l)));
                        u1_new[i](1, j, l) += rtskellam_cpp(
                            a_dis + b_dis * u1_new[i](1, j, l),
                            a_dis + b_dis * u1_new[i](1, j, l),
                            str1,
                            eng,
                            tempLB,
                            u1[i](0, j, l) - u1_new[i](1, j, l) + u1[i](1, j, l) - Ainc(j, l) - Pinc(j, l)
                        );
                        Einc(j, l) = u1_new[i](1, j, l) - u1[i](1, j, l) + Ainc(j, l) + Pinc(j, l);
                        if(Einc(j, l) < 0) Rprintf("Einc = %d j = %d l = %d\n", Einc(j, l), j, l);
                        if(u1_new[i](1, j, l) < 0) Rprintf("E = %d j = %d l = %d\n", u1_new[i](1, j, l), j, l);
                    }
                }
                
                // S given later
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        u1_new[i](0, j, l) = u1[i](0, j, l) - Einc(j, l);
                        if(u1_new[i](0, j, l) < 0) Rprintf("S = %d j = %d l = %d\n", u1_new[i](0, j, l), j, l);
                    }
                }
            } else {
                // DH incidence
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        DHinc(j, l) = u1_new[i](11, j, l) - u1[i](11, j, l);
                        if(DHinc(j, l) < 0) Rprintf("DHinc = %d j = %d l = %d\n", DHinc(j, l), j, l);
                    }
                }
                
                // DI incidence
                for(j = 0; j < nages; j++) {
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        DIinc(j, l) = u1_new[i](6, j, l) - u1[i](6, j, l);
                        if(DIinc(j, l) < 0) Rprintf("DIinc = %d j = %d l = %d\n", DIinc(j, l), j, l);
                    }
                }
            }
            
            // generate data in correct format for observation error weights
            Dtempinc.zeros();
            for(l = 0; l < u1_moves.n_rows; l++) {
                for(j = 0; j < nages; j++) {
                    Dtempinc((arma::uword) j * nlads + u1_moves(l, 0) - 1) += DIinc(j, l);
                    Dtempinc((arma::uword) nlads * nages + j * nlads + u1_moves(l, 0) - 1) += DHinc(j, l);
                }
            }
            
            // calculate log observation error weights
            for(l = 0; l < Dtempinc.n_elem; l++) {
                tempdens(l) = ldtskellam_cpp(
                    obsInc(l) - Dtempinc(l),
                    a1 + b * Dtempinc(l),
                    a2 + b * Dtempinc(l),
                    str1,
                    -Dtempinc(l),
                    obsInc(l),
                    0
                );
            }
            weights(i) += log_sum_exp(tempdens, 0);
            
            // advance seed
            seeds((arma::uword) omp_get_thread_num()) = eng();
        }
        
        // calculate log-likelihood contribution
        ll += log_sum_exp(weights, 1);
        
        // if zero likelihood then return
        if(!arma::is_finite(ll)) {
            if(saveAll == 0) {
                return List::create(Named("ll") = ll);
            } else {
                return List::create(Named("ll") = ll, _["particles"] = out);
            }
        }
        
        // normalise weights
        wnorm = log_sum_exp(weights, 0);
        weights = exp(weights - wnorm);
        weights = weights / sum(weights);
        
        // resample
        for(i = 0; i < npart; i++) {
            l = (arma::uword) rmultinom_cpp(weights, engSerial);
            u1[i] = u1_new[l];
        }
        
        // copy in order to pass by reference
        for(i = 0; i < npart; i++) {
            u1_new[i] = u1[i];
        }
        
        // save particles if necessary
        if(saveAll != 0) {
            if(saveAll == 1) {
                for(i = 0; i < npart; i++) {
                    // extract just counts for DI and DH
                    u1_night_reduced.zeros();
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        for(j = 0; j < nages; j++) {
                            u1_night_reduced(0, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](6, j, l);
                            u1_night_reduced(1, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](11, j, l);
                        }
                    }
                    out[i + npart * (t + 1)] = u1_night_reduced;
                }
            } else {
                for(i = 0; i < npart; i++) {
                    u1_night_full.zeros();
                    for(l = 0; l < u1_moves.n_rows; l++) {
                        for(j = 0; j < nages; j++) {
                            for(k = 0; k < nclasses; k++) {
                                u1_night_full(k, j, (arma::uword) u1_moves(l, 0) - 1) += u1[i](k, j, l);
                            }
                        }
                    }
                    out[i + npart * (t + 1)] = u1_night_full;
                }
            }
        }
        
        //calculate block run time
        timer.step("");
        NumericVector res(timer);
        
        Rprintf("t = %d / %d time = %.2f secs \n", t + 1, ndays, (res[timer_cnt] / 1e9) - prev_time);
        
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
