#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp/Benchmark/Timer.h>

//function to draw from a multivariate normal
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat cholSigma)
{
    int ncols = cholSigma.n_cols;
    arma::mat Y = arma::randn(n, ncols);
    return arma::repmat(mu, 1, n).t() + Y * cholSigma;
}

// function to compute Cholesky decomposition and to deal with
// cases where it fails
arma::mat cholArma(arma::mat sigma, double *scale)
{
    bool success = false;
    int j = 0;
    arma::mat cholSigma (sigma.n_rows, sigma.n_cols);
    while(success == false && j < 20) {
        success = arma::chol(cholSigma, sigma);
        if(success == false) {
            sigma += arma::eye(sigma.n_rows, sigma.n_cols) * 1e-6;
        }
        j++;
    }
//    if(j > 1) {
//        Rprintf("chol %d\n", j);
//    }
    // return scale
    *scale = 1e-6 * (j - 1);
    if(success == false) {
        // // write to file as check
        // FILE *filemcmc;
        // char outfile [64] = "output_chol.txt";
        // //append current iteration to runinfo.txt
        // filemcmc = fopen(outfile, "w");
        // for(j = 0; j < sigma.n_rows; j++) {
        //     for(k = 0; k < sigma.n_cols; k++) {
        //         fprintf(filemcmc, "%f ", sigma(j, k));
        //     }
        //     fprintf(filemcmc, "\n");
        // }
        // fclose(filemcmc);
        Rcpp::stop("Error in Cholesky decomposition");
    }
    return cholSigma;
}


//function to calculate variance-covariance matrix from posterior samples
void calcPost(int i, int npars, 
                  arma::vec *tempmn, arma::mat *meanmat, 
                  arma::mat *meanmat1, arma::mat posterior, 
                  arma::mat *propcov)
{
    int j, k, m;

    //calculates variance-covariance matrix
    //first: means
    for(j = 0; j < npars; j++) {
        (*tempmn)[j] = 0;
        for(k = 0; k <= i; k++) {
            (*tempmn)[j] += posterior(k, j);
        }
        (*tempmn)[j] = (*tempmn)[j] / ((double) i + 1);
    }
    //matrix of product of means
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*meanmat)(j, k) = (*tempmn)[j] * (*tempmn)[k];
        }
    }
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*propcov)(j, k) = 0.0;
        }
    }
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            for(m = 0; m <= i; m++) {
                (*propcov)(j, k) += posterior(m, j) * posterior(m, k);
            }
            (*propcov)(j, k) -= (i + 1) * (*meanmat)(j, k);
            (*propcov)(j, k) = (*propcov)(j, k) / ((double) i);
        }
    }
//        Rprintf("Initial cov:\n");
//        for(k = 0; k < npars; k++)
//        {
//            for(j = 0; j < npars; j++) Rprintf("%f ", (*propcov)(j, k));
//            Rprintf("\n");
//        }
    return;
}

//function to update variance-covariance matrix based on current samples
void adaptUpdate(int i, int npars, 
                  arma::vec *tempmn, arma::mat *meanmat, 
                  arma::mat *meanmat1, arma::vec posterior, 
                  arma::mat *propcov)
{
    int j, k;

    //recursively update mean and covariance matrix
    for(j = 0; j < npars; j++) {
        (*tempmn)[j] = ((*tempmn)[j] * i + posterior[j]) / ((double) i + 1);
    }
    //new matrix of product of means
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*meanmat1)(j, k) = (*tempmn)[j] * (*tempmn)[k];
        }
    }
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*propcov)(j, k) = (((double) (i - 1)) / ((double) i)) * (*propcov)(j, k) + (1.0 / ((double) i)) * (i * (*meanmat)(j, k) - (i + 1) * (*meanmat1)(j, k) + posterior[j] * posterior[k]);
        }
    }
    for(j = 0; j < npars; j++) {
        for(k = 0; k < npars; k++) {
            (*meanmat)(j, k) = (*meanmat1)(j, k);
        }
    }
    return;
}

// function to update transition probabilities
void updateTransProbs (arma::mat *pXsum, arma::imat *nXsum, arma::vec *pars, int Npop, double pini) {

    double pH = (*pars)(0);
    double pHD = (*pars)(1);
    double pI1 = (*pars)(2);
    double pI1D = (*pars)(3);
    double pI1H = (*pars)(4);
    double pI2 = (*pars)(5);
    double pP = (*pars)(6);
    double pE = (*pars)(7);
    double pEP = (*pars)(8);
    double beta = (*pars)(9);
    double betaA = (*pars)(10);
    
    double TP = 1.0 / (-log(1.0 - pP));
    double TI1 = 1.0 / (-log(1.0 - pI1));
    double TI2 = 1.0 / (-log(1.0 - pI2));
    double TA = TP + TI1 + TI2;
    double pA = 1.0 - exp(-1.0 / TA);
    
    double norm, norm1;
                
    // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
    //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    
    // E to E
    (*pXsum)(1, 1) = log(1.0 - pE);
    // E to A
    (*pXsum)(1, 2) = log(1.0 - pEP) + log(pE);
    // E to P
    (*pXsum)(1, 4) = log(pEP) + log(pE);

    // A to A
    (*pXsum)(2, 2) = log(1.0 - pA);
    // A to RA
    (*pXsum)(2, 3) = log(pA);

    // RA to RA
    (*pXsum)(3, 3) = 0.0;

    // P to P
    (*pXsum)(4, 4) = log(1.0 - pP);
    // P to I1
    (*pXsum)(4, 5) = log(pP);

    // I1 to I1
    (*pXsum)(5, 5) = log(1.0 - pI1);
    // I1 to I2
    (*pXsum)(5, 6) = log(1.0 - pI1H - pI1D) + log(pI1);
    // I1 to DI
    (*pXsum)(5, 8) = log(pI1D) + log(pI1);
    // I1 to H
    (*pXsum)(5, 9) = log(pI1H) + log(pI1);

    // I2 to I2
    (*pXsum)(6, 6) = log(1.0 - pI2);
    // I2 to RI
    (*pXsum)(6, 7) = log(pI2);

    // RI to RI
    (*pXsum)(7, 7) = 0.0;

    // DI to DI
    (*pXsum)(8, 8) = 0.0;

    // H to H
    (*pXsum)(9, 9) = log(1.0 - pH);
    // H to RH
    (*pXsum)(9, 10) = log(1.0 - pHD) + log(pH);
    // H to DH
    (*pXsum)(9, 11) = log(pHD) + log(pH);

    // RH to RH
    (*pXsum)(10, 10) = 0.0;

    // DH to DH
    (*pXsum)(11, 11) = 0.0;

    // S to S
    (*pXsum)(0, 0) = (-beta * ((*nXsum)(0, 4) + (*nXsum)(0, 5) + (*nXsum)(0, 6)) - betaA * beta * (*nXsum)(0, 2)) / ((double) Npop);
    
    // S to E
    norm1 = ((*pXsum)(0, 0) < 0.0 ? 0.0:(*pXsum)(0, 0));
    norm1 = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 0) - norm1));
    norm = (norm1 > log(pini) ? norm1:log(pini));
    (*pXsum)(0, 1) = norm + log(exp(norm1 - norm) + exp(log(pini) - norm) - exp(norm1 + log(pini) - norm));
    norm1 = ((*pXsum)(0, 1) < 0.0 ? 0.0:(*pXsum)(0, 1));
    (*pXsum)(0, 0) = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 1) - norm1));
    
//    // check sum
//    Rprintf("%f\n", exp((*pXsum)(0, 0)) + exp((*pXsum)(0, 1)));
    return;
}

// log-likelihood function
double loglikelihood (double beta, double betaA, arma::mat *pXsum, arma::imat *nXsum, arma::imat *nXcum, arma::ivec *t, double pini, int Npop, arma::ivec *D) {
                
    // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
    //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11

    double ll = (*nXsum)(0, 0) * log(1.0 - pini) + (*nXsum)(0, 1) * log(pini);
                
    // observation process
    ll += R::dunif((*D)(0), (*nXcum)(0, 8) + (*nXcum)(0, 11) - 0.1, (*nXcum)(0, 8) + (*nXcum)(0, 11) + 0.1, 1);
    
    double norm, norm1;
    for(int j = 1; j < (*t).n_elem; j++) {
        
        // S to S
        ll += ((*nXsum)(j - 1, 0) - ((*nXcum)(j, 1) - (*nXcum)(j - 1, 1))) * (*pXsum)(0, 0);
        // S to E
        ll += ((*nXcum)(j, 1) - (*nXcum)(j - 1, 1)) * (*pXsum)(0, 1);
        
        // E to E
        ll += ((*nXsum)(j - 1, 1) - ((*nXcum)(j, 2) - (*nXcum)(j - 1, 2)) - ((*nXcum)(j, 4) - (*nXcum)(j - 1, 4))) * (*pXsum)(1, 1);
        // E to A
        ll += ((*nXcum)(j, 2) - (*nXcum)(j - 1, 2)) * (*pXsum)(1, 2);
        // E to P
        ll += ((*nXcum)(j, 4) - (*nXcum)(j - 1, 4)) * (*pXsum)(1, 4);
        
        // A to A
        ll += ((*nXsum)(j - 1, 2) - ((*nXcum)(j, 3) - (*nXcum)(j - 1, 3))) * (*pXsum)(2, 2);
        // A to RA
        ll += ((*nXcum)(j, 3) - (*nXcum)(j - 1, 3)) * (*pXsum)(2, 3);
        
        // P to P
        ll += ((*nXsum)(j - 1, 4) - ((*nXcum)(j, 5) - (*nXcum)(j - 1, 5))) * (*pXsum)(4, 4);
        // P to I1
        ll += ((*nXcum)(j, 5) - (*nXcum)(j - 1, 5)) * (*pXsum)(4, 5);
        
        // I1 to I1
        ll += ((*nXsum)(j - 1, 5) - ((*nXcum)(j, 6) - (*nXcum)(j - 1, 6)) - ((*nXcum)(j, 8) - (*nXcum)(j - 1, 8)) - ((*nXcum)(j, 9) - (*nXcum)(j - 1, 9))) * (*pXsum)(5, 5);
        // I1 to I2
        ll += ((*nXcum)(j, 6) - (*nXcum)(j - 1, 6)) * (*pXsum)(5, 6);
        // I1 to DI
        ll += ((*nXcum)(j, 8) - (*nXcum)(j - 1, 8)) * (*pXsum)(5, 8);
        // I1 to H
        ll += ((*nXcum)(j, 9) - (*nXcum)(j - 1, 9)) * (*pXsum)(5, 9);
        
        // I2 to I2
        ll += ((*nXsum)(j - 1, 6) - ((*nXcum)(j, 7) - (*nXcum)(j - 1, 7))) * (*pXsum)(6, 6);
        // I2 to RI
        ll += ((*nXcum)(j, 7) - (*nXcum)(j - 1, 7)) * (*pXsum)(6, 7);
        
        // H to H
        ll += ((*nXsum)(j - 1, 9) - ((*nXcum)(j, 10) - (*nXcum)(j - 1, 10)) - ((*nXcum)(j, 11) - (*nXcum)(j - 1, 11))) * (*pXsum)(9, 9);
        // H to RH
        ll += ((*nXcum)(j, 10) - (*nXcum)(j - 1, 10)) * (*pXsum)(9, 10);
        // H to DH
        ll += ((*nXcum)(j, 11) - (*nXcum)(j - 1, 11)) * (*pXsum)(9, 11);
                
        // observation process
        ll += R::dunif((*D)(j), (*nXcum)(j, 8) + (*nXcum)(j, 11) - 0.1, (*nXcum)(j, 8) + (*nXcum)(j, 11) + 0.1, 1);
        
        // S to S
        (*pXsum)(0, 0) = (-beta * ((*nXsum)(j, 4) + (*nXsum)(j, 5) + (*nXsum)(j, 6)) - betaA * beta * (*nXsum)(j, 2)) / ((double) Npop);
        
        // S to E
        norm1 = ((*pXsum)(0, 0) < 0.0 ? 0.0:(*pXsum)(0, 0));
        norm1 = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 0) - norm1));
        norm = (norm1 > log(pini) ? norm1:log(pini));
        (*pXsum)(0, 1) = norm + log(exp(norm1 - norm) + exp(log(pini) - norm) - exp(norm1 + log(pini) - norm));
        norm1 = ((*pXsum)(0, 1) < 0.0 ? 0.0:(*pXsum)(0, 1));
        (*pXsum)(0, 0) = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 1) - norm1));
    }
    return ll;
}

// forwards filtering function
void forwardsFilter(int i, double beta, double betaA, arma::imat *Xind, arma::cube *pXind, arma::imat *nXsum, arma::imat *nXcum, arma::mat *pXsum, arma::ivec *t, arma::ivec *D, int Npop, double pini, arma::vec *tempProb, int print) {

    //set up variables
    int k, j, l;
    double norm, norm1;
                
    // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
    //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11

    // set initial state probabilities
    for(k = 0; k < 12; k++) {
        (*pXind)(i, 0, k) = R_NegInf;
    }
    (*pXind)(i, 0, 1) = log(pini);
    (*pXind)(i, 0, 0) = log(1.0 - pini);
    
    // print outputs if required
    if(print == 1) {
        Rprintf("Xind -1:\n");
        for(int r = 0; r < (*t).n_elem; r++) Rprintf("%d ", (*Xind)(i, r));
        Rprintf("\n");
        Rprintf("nXsum:\n");
        for(int r = 0; r < (*t).n_elem; r++) {
            for(int c = 0; c < 12; c++) {
                Rprintf("%d ", (*nXsum)(r, c));
            }
            Rprintf("\n");
        }
        
        Rprintf("nXcum:\n");
        for(int r = 0; r < (*t).n_elem; r++) {
            for(int c = 1; c < 12; c++) {
                Rprintf("%d ", (*nXcum)(r, c));
            }
            Rprintf("\n");
        }
    }
                
    // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
    //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    
    // loop over possible states
    for(k = 0; k < 12; k++) {
        
        if((*pXind)(i, 0, k) > R_NegInf) {
            
            // set observation process
            l = (*nXcum)(0, 8) + (k == 8 ? 1:0);
            l += (*nXcum)(0, 11) + (k == 11 ? 1:0);
            (*pXind)(i, 0, k) += R::dunif((*D)(0), l - 0.1, l + 0.1, 1);
                
            // if infective state then update transmission probs
            if(k == 4 || k == 5 || k == 6) {
                // update transmission probabilities
                // S to S
                (*pXsum)(0, 0) = (-beta * (1 + (*nXsum)(0, 4) + (*nXsum)(0, 5) + (*nXsum)(0, 6)) - betaA * beta * (*nXsum)(0, 2)) / ((double) Npop);
            } else {
                // if asymptomatic state then update transmission probs
                if(k == 2) {
                    // update transmission probabilities
                    // S to S
                    (*pXsum)(0, 0) = (-beta * ((*nXsum)(0, 4) + (*nXsum)(0, 5) + (*nXsum)(0, 6)) - betaA * beta * (1 + (*nXsum)(0, 2))) / ((double) Npop);
                } else {
                    // update transmission probabilities
                    // S to S
                    (*pXsum)(0, 0) = (-beta * ((*nXsum)(0, 4) + (*nXsum)(0, 5) + (*nXsum)(0, 6)) - betaA * beta * (*nXsum)(0, 2)) / ((double) Npop);
                }
            }
            // S to E
            norm1 = ((*pXsum)(0, 0) < 0.0 ? 0.0:(*pXsum)(0, 0));
            norm1 = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 0) - norm1));
            norm = (norm1 > log(pini) ? norm1:log(pini));
            (*pXsum)(0, 1) = norm + log(exp(norm1 - norm) + exp(log(pini) - norm) - exp(norm1 + log(pini) - norm));
            norm1 = ((*pXsum)(0, 1) < 0.0 ? 0.0:(*pXsum)(0, 1));
            (*pXsum)(0, 0) = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 1) - norm1));
            
            // now calculate other individual's transition probabilities
            
            // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
            //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
            
            // S to S
            (*pXind)(i, 0, k) += ((*nXsum)(0, 0) - ((*nXcum)(1, 1) - (*nXcum)(0, 1))) * (*pXsum)(0, 0);
            // S to E
            (*pXind)(i, 0, k) += ((*nXcum)(1, 1) - (*nXcum)(0, 1)) * (*pXsum)(0, 1);
            
            // E to E
            (*pXind)(i, 0, k) += ((*nXsum)(0, 1) - ((*nXcum)(1, 2) - (*nXcum)(0, 2)) - ((*nXcum)(1, 4) - (*nXcum)(0, 4))) * (*pXsum)(1, 1);
            // E to A
            (*pXind)(i, 0, k) += ((*nXcum)(1, 2) - (*nXcum)(0, 2)) * (*pXsum)(1, 2);
            // E to P
            (*pXind)(i, 0, k) += ((*nXcum)(1, 4) - (*nXcum)(0, 4)) * (*pXsum)(1, 4);
            
            // A to A
            (*pXind)(i, 0, k) += ((*nXsum)(0, 2) - ((*nXcum)(1, 3) - (*nXcum)(0, 3))) * (*pXsum)(2, 2);
            // A to RA
            (*pXind)(i, 0, k) += ((*nXcum)(1, 3) - (*nXcum)(0, 3)) * (*pXsum)(2, 3);
            
            // RA to RA
            (*pXind)(i, 0, k) += (*nXsum)(0, 3) * (*pXsum)(3, 3);
            
            // P to P
            (*pXind)(i, 0, k) += ((*nXsum)(0, 4) - ((*nXcum)(1, 5) - (*nXcum)(0, 5))) * (*pXsum)(4, 4);
            // P to I1
            (*pXind)(i, 0, k) += ((*nXcum)(1, 5) - (*nXcum)(0, 5)) * (*pXsum)(4, 5);
            
            // I1 to I1
            (*pXind)(i, 0, k) += ((*nXsum)(0, 5) - ((*nXcum)(1, 6) - (*nXcum)(0, 6)) - ((*nXcum)(1, 8) - (*nXcum)(0, 8)) - ((*nXcum)(1, 9) - (*nXcum)(0, 9))) * (*pXsum)(5, 5);
            // I1 to I2
            (*pXind)(i, 0, k) += ((*nXcum)(1, 6) - (*nXcum)(0, 6)) * (*pXsum)(5, 6);
            // I1 to DI
            (*pXind)(i, 0, k) += ((*nXcum)(1, 8) - (*nXcum)(0, 8)) * (*pXsum)(5, 8);
            // I1 to H
            (*pXind)(i, 0, k) += ((*nXcum)(1, 9) - (*nXcum)(0, 9)) * (*pXsum)(5, 9);
            
            // I2 to I2
            (*pXind)(i, 0, k) += ((*nXsum)(0, 6) - ((*nXcum)(1, 7) - (*nXcum)(0, 7))) * (*pXsum)(6, 6);
            // I2 to RI
            (*pXind)(i, 0, k) += ((*nXcum)(1, 7) - (*nXcum)(0, 7)) * (*pXsum)(6, 7);
            
            // RI to RI
            (*pXind)(i, 0, k) += (*nXsum)(0, 7) * (*pXsum)(7, 7);
            
            // H to H
            (*pXind)(i, 0, k) += ((*nXsum)(0, 9) - ((*nXcum)(1, 10) - (*nXcum)(0, 10)) - ((*nXcum)(1, 11) - (*nXcum)(0, 11))) * (*pXsum)(9, 9);
            // H to RH
            (*pXind)(i, 0, k) += ((*nXcum)(1, 10) - (*nXcum)(0, 10)) * (*pXsum)(9, 10);
            // H to DH
            (*pXind)(i, 0, k) += ((*nXcum)(1, 11) - (*nXcum)(0, 11)) * (*pXsum)(9, 11);
            
            // RH to RH
            (*pXind)(i, 0, k) += (*nXsum)(0, 10) * (*pXsum)(10, 10);
            
            // DH to DH
            (*pXind)(i, 0, k) += (*nXsum)(0, 11) * (*pXsum)(11, 11);
            
            if(std::isnan((*pXind)(i, 0, k))) Rcpp::stop("filter nan\n");
        }
    }
    
    // normalise probabilities on the log-scale
    norm = (*pXind)(i, 0, 0);
    for(k = 1; k < 12; k++) {
        norm = (norm > (*pXind)(i, 0, k) ? norm:(*pXind)(i, 0, k));
    }
    norm1 = 0.0;
    for(k = 0; k < 12; k++) {
        norm1 += exp((*pXind)(i, 0, k) - norm);
    }
    norm += log(norm1);
    for(k = 0; k < 12; k++) {
        (*pXind)(i, 0, k) -= norm;
    }
    
    // print outputs if required
    if(print == 1) {
        Rprintf("pXind(%d, 0): ", i);
        for(int c = 0; c < 12; c++) {
            Rprintf("%f ", exp((*pXind)(i, 0, c)));
        }
        Rprintf("\n");
    }
    
    // loop over other time points
    for(j = 1; j < ((*t).n_elem - 1); j++) {
        
        // loop over possible states to move to
        for(l = 0; l < 12; l++) {
                
            // loop over possible states to move from
            for(k = 0; k < 12; k++) {
            
                if((*pXind)(i, j - 1, k) > R_NegInf) {
                
                    if(l == 0 || l == 1) {
                        // update transmission probs
                        if(k == 4 || k == 5 || k == 6) {
                            // update transmission probabilities
                            // S to S
                            (*pXsum)(0, 0) = (-beta * (1 + (*nXsum)(j - 1, 4) + (*nXsum)(j - 1, 5) + (*nXsum)(j - 1, 6)) - betaA * beta * (*nXsum)(j - 1, 2)) / ((double) Npop);
                        } else {
                            // if asymptomatic state then update transmission probs
                            if(k == 2) {
                                // update transmission probabilities
                                // S to S
                                (*pXsum)(0, 0) = (-beta * ((*nXsum)(j - 1, 4) + (*nXsum)(j - 1, 5) + (*nXsum)(j - 1, 6)) - betaA * beta * (1 + (*nXsum)(j - 1, 2))) / ((double) Npop);
                            } else {
                                // update transmission probabilities
                                // S to S
                                (*pXsum)(0, 0) = (-beta * ((*nXsum)(j - 1, 4) + (*nXsum)(j - 1, 5) + (*nXsum)(j - 1, 6)) - betaA * beta * (*nXsum)(j - 1, 2)) / ((double) Npop);
                            }
                        }
                        // S to E
                        norm1 = ((*pXsum)(0, 0) < 0.0 ? 0.0:(*pXsum)(0, 0));
                        norm1 = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 0) - norm1));
                        norm = (norm1 > log(pini) ? norm1:log(pini));
                        (*pXsum)(0, 1) = norm + log(exp(norm1 - norm) + exp(log(pini) - norm) - exp(norm1 + log(pini) - norm));
                        norm1 = ((*pXsum)(0, 1) < 0.0 ? 0.0:(*pXsum)(0, 1));
                        (*pXsum)(0, 0) = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 1) - norm1));
                    }
                    
                    (*tempProb)(k) = (*pXind)(i, j - 1, k) + (*pXsum)(k, l);
                } else {
                    (*tempProb)(k) = R_NegInf;
                }
            }
            // sum over states to move from
            norm = (*tempProb)(0);
            for(k = 1; k < 12; k++) {
                norm = (norm > (*tempProb)(k) ? norm:(*tempProb)(k));
            }
            if(norm > R_NegInf) {
                (*pXind)(i, j, l) = 0.0;
                for(k = 0; k < 12; k++) {
                    (*pXind)(i, j, l) += exp((*tempProb)(k) - norm);
                }
                (*pXind)(i, j, l) = norm + log((*pXind)(i, j, l));
            } else {
                (*pXind)(i, j, l) = R_NegInf;
            }
        }
        
        // check
        norm = 0.0;
        for(k = 0; k < 12; k++) {
            norm += exp((*pXind)(i, j, k));
        }
        if(fabs(norm - 1.0) > 1e-5) {
            Rprintf("norm = %f\n", norm);
            Rcpp::stop("norm mismatch\n");
        }
        
        // loop over possible states to move to
        for(k = 0; k < 12; k++) {
            
            if((*pXind)(i, j, k) > R_NegInf) {
            
                // set observation process
                l = (*nXcum)(j, 8) + (k == 8 ? 1:0);
                l += (*nXcum)(j, 11) + (k == 11 ? 1:0);
                (*pXind)(i, j, k) += R::dunif((*D)(j), l - 0.1, l + 0.1, 1);
                    
                // if infective state then update transmission probs
                if(k == 4 || k == 5 || k == 6) {
                    // update transmission probabilities
                    // S to S
                    (*pXsum)(0, 0) = (-beta * (1 + (*nXsum)(j, 4) + (*nXsum)(j, 5) + (*nXsum)(j, 6)) - betaA * beta * (*nXsum)(j, 2)) / ((double) Npop);
                } else {
                    if(k == 2) {
                        // update transmission probabilities
                        // S to S
                        (*pXsum)(0, 0) = (-beta * ((*nXsum)(j, 4) + (*nXsum)(j, 5) + (*nXsum)(j, 6)) - betaA * beta * (1 + (*nXsum)(j, 2))) / ((double) Npop);
                    } else {
                        // update transmission probabilities
                        // S to S
                        (*pXsum)(0, 0) = (-beta * ((*nXsum)(j, 4) + (*nXsum)(j, 5) + (*nXsum)(j, 6)) - betaA * beta * (*nXsum)(j, 2)) / ((double) Npop);
                    }
                }
                // S to E
                norm1 = ((*pXsum)(0, 0) < 0.0 ? 0.0:(*pXsum)(0, 0));
                norm1 = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 0) - norm1));
                norm = (norm1 > log(pini) ? norm1:log(pini));
                (*pXsum)(0, 1) = norm + log(exp(norm1 - norm) + exp(log(pini) - norm) - exp(norm1 + log(pini) - norm));
                norm1 = ((*pXsum)(0, 1) < 0.0 ? 0.0:(*pXsum)(0, 1));
                (*pXsum)(0, 0) = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 1) - norm1));
                
                // now calculate other individual's transition probabilities
                
                // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
                //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
                
                // S to S
                (*pXind)(i, j, k) += ((*nXsum)(j, 0) - ((*nXcum)(j + 1, 1) - (*nXcum)(j, 1))) * (*pXsum)(0, 0);
                // S to E
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 1) - (*nXcum)(j, 1)) * (*pXsum)(0, 1);
                
                // E to E
                (*pXind)(i, j, k) += ((*nXsum)(j, 1) - ((*nXcum)(j + 1, 2) - (*nXcum)(j, 2)) - ((*nXcum)(j + 1, 4) - (*nXcum)(j, 4))) * (*pXsum)(1, 1);
                // E to A
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 2) - (*nXcum)(j, 2)) * (*pXsum)(1, 2);
                // E to P
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 4) - (*nXcum)(j, 4)) * (*pXsum)(1, 4);
                
                // A to A
                (*pXind)(i, j, k) += ((*nXsum)(j, 2) - ((*nXcum)(j + 1, 3) - (*nXcum)(j, 3))) * (*pXsum)(2, 2);
                // A to RA
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 3) - (*nXcum)(j, 3)) * (*pXsum)(2, 3);
                
                // RA to RA
                (*pXind)(i, j, k) += (*nXsum)(j, 3) * (*pXsum)(3, 3);
                
                // P to P
                (*pXind)(i, j, k) += ((*nXsum)(j, 4) - ((*nXcum)(j + 1, 5) - (*nXcum)(j, 5))) * (*pXsum)(4, 4);
                // P to I1
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 5) - (*nXcum)(j, 5)) * (*pXsum)(4, 5);
                
                // I1 to I1
                (*pXind)(i, j, k) += ((*nXsum)(j, 5) - ((*nXcum)(j + 1, 6) - (*nXcum)(j, 6)) - ((*nXcum)(j + 1, 8) - (*nXcum)(j, 8)) - ((*nXcum)(j + 1, 9) - (*nXcum)(j, 9))) * (*pXsum)(5, 5);
                // I1 to I2
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 6) - (*nXcum)(j, 6)) * (*pXsum)(5, 6);
                // I1 to DI
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 8) - (*nXcum)(j, 8)) * (*pXsum)(5, 8);
                // I1 to H
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 9) - (*nXcum)(j, 9)) * (*pXsum)(5, 9);
                
                // I2 to I2
                (*pXind)(i, j, k) += ((*nXsum)(j, 6) - ((*nXcum)(j + 1, 7) - (*nXcum)(j, 7))) * (*pXsum)(6, 6);
                // I2 to RI
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 7) - (*nXcum)(j, 7)) * (*pXsum)(6, 7);
                
                // RI to RI
                (*pXind)(i, j, k) += (*nXsum)(j, 7) * (*pXsum)(7, 7);
                
                // H to H
                (*pXind)(i, j, k) += ((*nXsum)(j, 9) - ((*nXcum)(j + 1, 10) - (*nXcum)(j, 10)) - ((*nXcum)(j + 1, 11) - (*nXcum)(j, 11))) * (*pXsum)(9, 9);
                // H to RH
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 10) - (*nXcum)(j, 10)) * (*pXsum)(9, 10);
                // H to DH
                (*pXind)(i, j, k) += ((*nXcum)(j + 1, 11) - (*nXcum)(j, 11)) * (*pXsum)(9, 11);
                
                // RH to RH
                (*pXind)(i, j, k) += (*nXsum)(j, 10) * (*pXsum)(10, 10);
                
                // DH to DH
                (*pXind)(i, j, k) += (*nXsum)(j, 11) * (*pXsum)(11, 11);
            }
        }
        
        // normalise probabilities on the log-scale
        norm = (*pXind)(i, j, 0);
        for(k = 1; k < 12; k++) {
            norm = (norm > (*pXind)(i, j, k) ? norm:(*pXind)(i, j, k));
        }
        norm1 = 0.0;
        for(k = 0; k < 12; k++) {
            norm1 += exp((*pXind)(i, j, k) - norm);
        }
        norm += log(norm1);
        for(k = 0; k < 12; k++) {
            (*pXind)(i, j, k) -= norm;
        }
        
        // print outputs if required
        if(print == 1) {
            Rprintf("pXind(%d, %d): ", i, j);
            for(int c = 0; c < 12; c++) {
                Rprintf("%f ", exp((*pXind)(i, j, c)));
            }
            Rprintf("\n");
        }
    }
    
    // update final time point
    
    // loop over possible states to move to
    for(l = 0; l < 12; l++) {
            
        // loop over possible states to move from
        for(k = 0; k < 12; k++) {
        
            if((*pXind)(i, j - 1, k) > R_NegInf) {
            
                if(l == 0 || l == 1) {
                    // update transmission probs
                    if(k == 4 || k == 5 || k == 6) {
                        // update transmission probabilities
                        // S to S
                        (*pXsum)(0, 0) = (-beta * (1 + (*nXsum)(j - 1, 4) + (*nXsum)(j - 1, 5) + (*nXsum)(j - 1, 6)) - betaA * beta * (*nXsum)(j - 1, 2)) / ((double) Npop);
                    } else {
                        // if asymptomatic state then update transmission probs
                        if(k == 2) {
                            // update transmission probabilities
                            // S to S
                            (*pXsum)(0, 0) = (-beta * ((*nXsum)(j - 1, 4) + (*nXsum)(j - 1, 5) + (*nXsum)(j - 1, 6)) - betaA * beta * (1 + (*nXsum)(j - 1, 2))) / ((double) Npop);
                        } else {
                            // update transmission probabilities
                            // S to S
                            (*pXsum)(0, 0) = (-beta * ((*nXsum)(j - 1, 4) + (*nXsum)(j - 1, 5) + (*nXsum)(j - 1, 6)) - betaA * beta * (*nXsum)(j - 1, 2)) / ((double) Npop);
                        }
                    }
                    // S to E
                    norm1 = ((*pXsum)(0, 0) < 0.0 ? 0.0:(*pXsum)(0, 0));
                    norm1 = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 0) - norm1));
                    norm = (norm1 > log(pini) ? norm1:log(pini));
                    (*pXsum)(0, 1) = norm + log(exp(norm1 - norm) + exp(log(pini) - norm) - exp(norm1 + log(pini) - norm));
                    norm1 = ((*pXsum)(0, 1) < 0.0 ? 0.0:(*pXsum)(0, 1));
                    (*pXsum)(0, 0) = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 1) - norm1));
                }
                
                (*tempProb)(k) = (*pXind)(i, j - 1, k) + (*pXsum)(k, l);
            } else {
                (*tempProb)(k) = R_NegInf;
            }
        }
        // sum over states to move from
        norm = (*tempProb)(0);
        for(k = 1; k < 12; k++) {
            norm = (norm > (*tempProb)(k) ? norm:(*tempProb)(k));
        }
        if(norm > R_NegInf) {
            (*pXind)(i, j, l) = 0.0;
            for(k = 0; k < 12; k++) {
                (*pXind)(i, j, l) += exp((*tempProb)(k) - norm);
            }
            (*pXind)(i, j, l) = norm + log((*pXind)(i, j, l));
        } else {
            (*pXind)(i, j, l) = R_NegInf;
        }
    }
        
    // check
    norm = 0.0;
    for(k = 0; k < 12; k++) {
        norm += exp((*pXind)(i, j, k));
    }
    if(fabs(norm - 1.0) > 1e-5) {
        Rprintf("norm = %f\n", norm);
        Rcpp::stop("norm mismatch1\n");
    }
    
    // loop over possible states to move to
    for(k = 0; k < 12; k++) {
        if((*pXind)(i, j, k) > R_NegInf) {
            // set observation process
            l = (*nXcum)(j, 8) + (k == 8 ? 1:0);
            l += (*nXcum)(j, 11) + (k == 11 ? 1:0);
            (*pXind)(i, j, k) += R::dunif((*D)(j), l - 0.1, l + 0.1, 1);
        }
    }  
    
    // normalise probabilities on the log-scale
    norm = (*pXind)(i, j, 0);
    for(k = 1; k < 12; k++) {
        norm = (norm > (*pXind)(i, j, k) ? norm:(*pXind)(i, j, k));
    }
    norm1 = 0.0;
    for(k = 0; k < 12; k++) {
        norm1 += exp((*pXind)(i, j, k) - norm);
    }
    norm += log(norm1);
    for(k = 0; k < 12; k++) {
        (*pXind)(i, j, k) -= norm;
    }
    
    // print outputs if required
    if(print == 1) {
        Rprintf("pXind(%d, %d): ", i, j);
        for(int c = 0; c < 12; c++) {
            Rprintf("%f ", exp((*pXind)(i, j, c)));
        }
        Rprintf("\n");
    }
    return;
}

// backwards sampling step
void backwardsSampler (int i, double beta, double betaA, arma::cube *pXind, arma::ivec *t, arma::imat *Xind, arma::mat *pXsum, arma::imat *nXsum, arma::vec *pXback, int Npop, double pini, int print) {
    
    double u = R::runif(0.0, 1.0);
    int j = 0, k = 0;
    double norm1;
    double norm = exp((*pXind)(i, (*t).n_elem - 1, k));
    while(u > norm) {
        k++;
        norm += exp((*pXind)(i, (*t).n_elem - 1, k));
    }
    // set sample for individual i
    (*Xind)(i, (*t).n_elem - 1) = k;
    
    // print outputs if required
    if(print == 1) {
        Rprintf("\npXback(%d): ", (*t).n_elem - 1);
        for(int c = 0; c < 12; c++) {
            Rprintf("%f ", exp((*pXind)(i, (*t).n_elem - 1, c)));
        }
        Rprintf("\n");
    }
    
    // run through rest of time points
    for(j = ((*t).n_elem - 2); j >= 0; j--) {
        
        // loop over possible states at time j
        for(k = 0; k < 12; k++) {
        
            if((*Xind)(i, j + 1) == 0 || (*Xind)(i, j + 1) == 1) {
                // update transmission probs
                if(k == 4 || k == 5 || k == 6) {
                    // update transmission probabilities
                    // S to S
                    (*pXsum)(0, 0) = (-beta * (1 + (*nXsum)(j, 4) + (*nXsum)(j, 5) + (*nXsum)(j, 6)) - betaA * beta * (*nXsum)(j, 2)) / ((double) Npop);
                } else {
                    // if asymptomatic state then update transmission probs
                    if(k == 2) {
                        // update transmission probabilities
                        // S to S
                        (*pXsum)(0, 0) = (-beta * ((*nXsum)(j, 4) + (*nXsum)(j, 5) + (*nXsum)(j, 6)) - betaA * beta * (1 + (*nXsum)(j, 2))) / ((double) Npop);
                    } else {
                        // update transmission probabilities
                        // S to S
                        (*pXsum)(0, 0) = (-beta * ((*nXsum)(j, 4) + (*nXsum)(j, 5) + (*nXsum)(j, 6)) - betaA * beta * (*nXsum)(j, 2)) / ((double) Npop);
                    }
                }
                // S to E
                norm1 = ((*pXsum)(0, 0) < 0.0 ? 0.0:(*pXsum)(0, 0));
                norm1 = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 0) - norm1));
                norm = (norm1 > log(pini) ? norm1:log(pini));
                (*pXsum)(0, 1) = norm + log(exp(norm1 - norm) + exp(log(pini) - norm) - exp(norm1 + log(pini) - norm));
                norm1 = ((*pXsum)(0, 1) < 0.0 ? 0.0:(*pXsum)(0, 1));
                (*pXsum)(0, 0) = norm1 + log(exp(-norm1) - exp((*pXsum)(0, 1) - norm1));
            }
            
            // calculate unnormalised event probabilities
            (*pXback)(k) = (*pXsum)(k, (*Xind)(i, j + 1)) + (*pXind)(i, j, k);
        }
        
        // normalise
        norm = (*pXback)(0);
        for(k = 1; k < 12; k++) {
            norm = (norm > (*pXback)(k) ? norm:(*pXback)(k));
        } 
        if(norm > R_NegInf) {
            norm1 = 0.0;
            for(k = 0; k < 12; k++) {
                norm1 += exp((*pXback)(k) - norm);
            }
            norm += log(norm1);
            for(k = 0; k < 12; k++) {
                (*pXback)(k) -= norm;
            }
        } else {
            Rcpp::stop("Worrying norm1\n");
        }
        
        // print outputs if required
        if(print == 1) {
            Rprintf("pXback(%d): ", j);
            for(int c = 0; c < 12; c++) {
                Rprintf("%f ", exp((*pXback)(c)));
            }
            Rprintf("\n");
        }
        
        // sample state
        u = R::runif(0.0, 1.0);
        k = 0;
        norm = exp((*pXback)(k));
        while(u > norm) {
            k++;
            norm += exp((*pXback)(k));
        }
        // set sample for individual i
        (*Xind)(i, j) = k;
    }
    
    // print outputs if required
    if(print == 1) {
        Rprintf("Xind(%d):\n", i);
        for(j = 0; j < (*t).n_elem; j++) {
            Rprintf("%d ", (*Xind)(i, j));
        }
        Rprintf("\n");
    }
    return;
}

// [[Rcpp::export]]
arma::mat iFFBS(
    arma::ivec D_prime,
    arma::ivec t,
    int Npop,
    int Niter,
    double pini,
    arma::vec pars,
    int fixPars = 0,
    int outputCum = 0,
    int print = 0
) {
    // initialise counters and temporary vars
    int i, j, k, l, iter, Ntemp;
    double norm, norm1, u, psus;
    
    // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
    //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    
    // calculate how many infections (for initialisation)
    int Ninf = 0, Ninf_prev = 0;
    for(j = 0; j < t.n_elem; j++) {
        Ninf += D_prime(j);
    }
    int Nsus = Npop - Ninf;
    int Nposs = 10 * Ninf + 1;
    Nposs = (Nposs > Npop ? Npop:Nposs);
    if(print == 1) {
        Rprintf("Npop = %d Ninf = %d Nposs = %d\n", Npop, Ninf, Nposs);
    }
    
    // set up observed cumulative counts
    arma::ivec D(D_prime.n_elem);
    D(0) = D_prime(0);
    for(i = 1; i < D_prime.n_elem; i++) D(i) = D(i - 1) + D_prime(i);
    
    // print outputs if required
    if(print == 1) {
        Rprintf("D:\n");
        for(int c = 0; c < t.n_elem; c++) {
            Rprintf("%d ", D(c));
        }
        Rprintf("\n");
    }
    
    // calculate cumulative counts
    // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
    //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
    arma::imat nXcum;
    nXcum.zeros(t.n_elem, 12);
    
    // initialise number in states
    arma::imat nXsum;
    nXsum.zeros(t.n_elem, 12);
    for(j = 0; j < t.n_elem; j++) {
        nXsum(j, 0) = Npop;
    }
    
    // initialise transition probabilities
    arma::mat pXsum(12, 12);
    for(i = 0; i < 12; i++) {
        for(j = 0; j < 12; j++) {
            pXsum(i, j) = R_NegInf;
        }
    }
    
    // initialise filtering probabilities
    arma::cube pXind(Nposs, t.n_elem, 12);
    for(i = 0; i < Nposs; i++) {
        for(j = 0; j < t.n_elem; j++) {
            for(k = 0; k < 12; k++) {
                pXind(i, j, k) = R_NegInf;
            }
        }
    }
    
    // initialise backwards sampling probabilities
    arma::vec pXback;
    pXback.zeros(12);
    
    // state vectors
    arma::imat Xind;
    Xind.zeros(Nposs, t.n_elem);
    
    // set up temporary vector for filter
    arma::vec tempProb(12);
    
    // sample initial values from priors
    double acc = NAN, accprop, accrate;
    int setini = (pars.is_finite() ? 1:0);
    if(fixPars == 1 && setini != 1) {
        Rcpp::stop("Must set parameters if 'fixPars == 1'\n");
    }
    while(std::isnan(acc)) {
    
        // check for user interrupt
        R_CheckUserInterrupt();
        
        if(setini == 0) {
            pars(0) = R::runif(0.0, 1.0);
            pars(1) = R::runif(0.0, 1.0);
            pars(2) = R::runif(0.0, 1.0);
            pars(3) = R::runif(0.0, 1.0);
            pars(4) = R::runif(0.0, pars(3));
            pars(5) = R::runif(0.0, 1.0);
            pars(6) = R::runif(0.0, 1.0);
            pars(7) = R::runif(0.0, 1.0);
            pars(8) = R::runif(0.0, 1.0);
            pars(9) = R::rgamma(1.0, 1.0);
            pars(10) = R::runif(0.0, 1.0);
        }
        
        // initialise states so they are consistent
        // with the observed data
        Xind.zeros();
        k = 0;
        for(j = (t.n_elem - 1); j >= 0; j--) {
            for(i = 0; i < D_prime(j); i++) {
                if(R::runif(0.0, 1.0) < 0.5) {
                    Xind(k, j) = 8;
                    Xind(k, j - 1) = 5;
                    Xind(k, j - 2) = 4;
                    Xind(k, j - 3) = 1;
                } else {
                    Xind(k, j) = 11;
                    Xind(k, j - 1) = 9;
                    Xind(k, j - 2) = 5;
                    Xind(k, j - 3) = 4;
                    Xind(k, j - 4) = 1;
                } 
                k++;
            }
        }
        for(i = 0; i < Ninf; i++) {
            for(j = 1; j < t.n_elem; j++) {
                if(Xind(i, j) != Xind(i, j - 1) && Xind(i, j) == 0) {
                    Xind(i, j) = Xind(i, j - 1);
                }
            }
        }
        
        // calculate cumulative counts
        // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
        //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
        
        for(i = 0; i < Ninf; i++) {
            if(Xind(i, 0) > 0) {
                nXcum(0, Xind(i, 0))++;
            }
        }
        for(j = 1; j < t.n_elem; j++) {
            for(i = 0; i < 12; i++) {
                nXcum(j, i) = nXcum(j - 1, i);
            }
            for(i = 0; i < Ninf; i++) {
                if(Xind(i, j - 1) != Xind(i, j)) {
                    nXcum(j, Xind(i, j))++;
                }
            }
        }
        
        // initialise number in states
        for(i = 0; i < Ninf; i++) {
            for(j = 0; j < t.n_elem; j++) {
                nXsum(j, Xind(i, j))++;
            }
        }
        for(j = 0; j < t.n_elem; j++) {
            nXsum(j, 0) = Npop - nXcum(j, 1);
        }
        
        // print outputs if required
        if(print == 1) {
            Rprintf("nXsum:\n");
            for(int r = 0; r < t.n_elem; r++) {
                for(int c = 0; c < 12; c++) {
                    Rprintf("%d ", nXsum(r, c));
                }
                Rprintf("\n");
            }
        }
        
        // print outputs if required
        if(print == 1) {
            Rprintf("nXcum:\n");
            for(int r = 0; r < t.n_elem; r++) {
                for(int c = 1; c < 12; c++) {
                    Rprintf("%d ", nXcum(r, c));
                }
                Rprintf("\n");
            }
        }
        
        // update transition probabilities
        updateTransProbs(&pXsum, &nXsum, &pars, Npop, pini);
        
        // print outputs if required
        if(print == 1) {
            Rprintf("pXsum:\n");
            for(int r = 0; r < 12; r++) {
                for(int c = 0; c < 12; c++) {
                    Rprintf("%f ", pXsum(r, c));
                }
                Rprintf("\n");
            }
        }
        
        // calculate log-likelihood
        acc = loglikelihood(pars(9), pars(10), &pXsum, &nXsum, &nXcum, &t, pini, Npop, &D);
        
        // add priors
        acc += R::dexp(pars(9), 1.0, 1);
        
        // check for valid posterior
        if(setini == 1 && std::isnan(acc)) {
            Rcpp::stop("Can't initialise with these initial values\n");
        }
    }
    
    if(print == 1) {
        Rprintf("Xind:\n");
        for(int r = 0; r < Nposs; r++) {
            for(int c = 0; c < t.n_elem; c++) {
                Rprintf("%d ", Xind(r, c));
            }
            Rprintf("\n");
        }
    }
        
    // print outputs if required
    if(print == 1) {
        Rprintf("pXsum:\n");
        for(int r = 0; r < 12; r++) {
            for(int c = 0; c < 12; c++) {
                Rprintf("%f ", pXsum(r, c));
            }
            Rprintf("\n");
        }
    }
    
    // print outputs if required
    if(print == 1) {
        Rprintf("nXsum:\n");
        for(int r = 0; r < t.n_elem; r++) {
            for(int c = 0; c < 12; c++) {
                Rprintf("%d ", nXsum(r, c));
            }
            Rprintf("\n");
        }
    }
    
    // print outputs if required
    if(print == 1) {
        Rprintf("nXcum:\n");
        for(int r = 0; r < t.n_elem; r++) {
            for(int c = 1; c < 12; c++) {
                Rprintf("%d ", nXcum(r, c));
            }
            Rprintf("\n");
        }
    }
    
    // set up output vector
    int npars = 11;
    arma::mat output;
    if(fixPars == 0) {
        output.zeros(Niter, npars + t.n_elem * 12 + 1);
    } else {
        output.zeros(Niter, t.n_elem * 12 + 1);
    }
    
    // initialise chain and set up vector to hold proposals
    arma::vec parsProp(npars);
    parsProp.zeros();
    
    // set up adaptive proposal distribution
    double cholScale = 0.0;
    arma::mat propVar(npars, npars);
    arma::mat propVarIni(npars, npars);
    arma::mat propVarChol(npars, npars);
    arma::mat propVarIniChol(npars, npars);
    propVar.zeros(); propVarIni.zeros();
    for(j = 0; j < npars; j++) {
        propVarIni(j, j) = pow(0.1, 2.0) / ((double) npars);
        propVar(j, j) = pow(0.1, 2.0) / ((double) npars);
    }
    // calculate Cholesky decomposition for MVN sampling
    double adaptscale = pow(2.38, 2.0) / ((double) npars);
    propVarIniChol = cholArma(propVarIni, &cholScale);
    propVarChol = cholArma(propVar * adaptscale, &cholScale);
    arma::vec tempmn(npars);
    arma::mat meanmat(npars, npars);
    arma::mat meanmat1(npars, npars);
    
    //initialise timer
    Rcpp::Timer timer;
    int timer_cnt = 0;
    double prev_time = 0.0;
    
    // set up acceptance counters
    int nacc = 0, cumacc = 0;
    
    // checks
    arma::imat nXcum1;
    nXcum1.zeros(t.n_elem, 12);
    
    // initialise number in states
    arma::imat nXsum1;
    nXsum1.zeros(t.n_elem, 12);
    for(j = 0; j < t.n_elem; j++) {
        nXsum1(j, 0) = Npop;
    }
    
    // loop over iterations
    Rprintf("Starting MCMC...\n");
    for(iter = 0; iter < Niter; iter++) {
    
        // check for user interrupt
        R_CheckUserInterrupt();
        
        /////////////////////////////////////////////////////////////////
        ///////                     CHECKS                        ///////
        /////////////////////////////////////////////////////////////////
        
        // calculate cumulative counts
        // states: S, E, A, RA, P, I1, I2, RI, DI, H, RH, DH
        //         0, 1, 2, 3,  4, 5,  6,  7,  8,  9, 10, 11
        nXcum1.zeros();
        nXsum1.zeros();
        for(i = 0; i < Ninf; i++) {
            if(Xind(i, 0) > 0) {
                nXcum1(0, Xind(i, 0))++;
            }
        }
        for(j = 1; j < t.n_elem; j++) {
            for(i = 0; i < 12; i++) {
                nXcum1(j, i) = nXcum1(j - 1, i);
            }
            for(i = 0; i < Ninf; i++) {
                if(Xind(i, j - 1) != Xind(i, j)) {
                    nXcum1(j, Xind(i, j))++;
                }
            }
        }
        
        // initialise number in states
        for(i = 0; i < Ninf; i++) {
            for(j = 0; j < t.n_elem; j++) {
                nXsum1(j, Xind(i, j))++;
            }
        }
        for(j = 0; j < t.n_elem; j++) {
            nXsum1(j, 0) = Npop - nXcum1(j, 1);
        }
        for(i = 0; i < 12; i++) {
            for(j = 0; j < t.n_elem; j++) {
                if(nXsum(j, i) != nXsum1(j, i) || (nXcum(j, i) != nXcum1(j, i) && i > 0)) {
                    Rprintf("%d %d cum1 = %d cum = %d sum1 = %d sum = %d\n", j, i, nXcum1(j, i), nXcum(j, i), nXsum1(j, i), nXsum(j, i));
                    Rcpp::stop("Some error in counts\n");
                }
            }
        }
        
        /////////////////////////////////////////////////////////////////
        ///////                   PARAMETERS                      ///////
        /////////////////////////////////////////////////////////////////
        
        if(fixPars == 0) {
            // update transition probabilities
            updateTransProbs(&pXsum, &nXsum, &pars, Npop, pini);
            
            // propose new values
            if(R::runif(0.0, 1.0) < 0.05) {
                parsProp = arma::conv_to<arma::vec>::from(mvrnormArma(1, pars, propVarIniChol));
            } else {
                parsProp = arma::conv_to<arma::vec>::from(mvrnormArma(1, pars, propVarChol));
            }
            
            // check validity
            k = 0;
            for(j = 0; j < 9; j++) {
                k += (parsProp(j) > 0.0 && parsProp(j) < 1.0 ? 0:1);
            }
            k += ((parsProp(3) + parsProp(4)) < 1.0 ? 0:1);
            k += (parsProp(9) > 0.0 ? 0:1);
            k += (parsProp(10) > 0.0 && parsProp(10) < 1.0 ? 0:1);
            
            // if valid then do accept-reject step, else reject immediately
            if(k == 0) {  
                
                // calculate log-likelihood
                acc = loglikelihood(pars(9), pars(10), &pXsum, &nXsum, &nXcum, &t, pini, Npop, &D);
                
                // add priors
                acc += R::dexp(pars(9), 1.0, 1);
            
                // update transition probabilities
                updateTransProbs(&pXsum, &nXsum, &parsProp, Npop, pini);    
                
                // calculate log-likelihood
                accprop = loglikelihood (parsProp(9), parsProp(10), &pXsum, &nXsum, &nXcum, &t, pini, Npop, &D);
                
                // add priors
                accprop += R::dexp(parsProp(9), 1.0, 1);
                
                // calculate acceptance probability
                accrate = accprop - acc;
    //            Rprintf("accrate = %f\n", accrate);
                if(log(R::runif(0.0, 1.0)) < accrate) {
                    pars = parsProp;
                    nacc++;
                    cumacc++;
                    acc = accprop;
                }
            }
        }
        
        /////////////////////////////////////////////////////////////////
        ///////                    INFECTED                       ///////
        /////////////////////////////////////////////////////////////////
        
//        Rprintf("inf\n");
        
        // loop over individuals
        Ninf_prev = Ninf;
        i = 0;
        while(i < Ninf) {
        
            // update transition probabilities
            updateTransProbs(&pXsum, &nXsum, &pars, Npop, pini);  
            
            // remove individual i from states
            // and cumulative states to get correct
            // powers in subsequent updates
            nXsum(0, Xind(i, 0))--;
            for(k = 0; k < t.n_elem; k++) {
                nXcum(k, Xind(i, 0))--;
            }
            for(k = 1; k < t.n_elem; k++) {
                if(Xind(i, k - 1) != Xind(i, k)) {
                    for(j = k; j < t.n_elem; j++) {
                        nXcum(j, Xind(i, k))--;
                    }
                }   
                nXsum(k, Xind(i, k))--;
            }
            
            // forwards filter
            forwardsFilter(i, pars(9), pars(10), &Xind, &pXind, &nXsum, &nXcum, &pXsum, &t, &D, Npop, pini, &tempProb, print);
            
            // backwards sampler
            backwardsSampler (i, pars(9), pars(10), &pXind, &t, &Xind, &pXsum, &nXsum, &pXback, Npop, pini, print);
            
            // add individual i from states
            // and cumulative states to get correct
            // powers in subsequent updates
            nXsum(0, Xind(i, 0))++;
            for(k = 0; k < t.n_elem; k++) {
                nXcum(k, Xind(i, 0))++;
            }
            for(k = 1; k < t.n_elem; k++) {
                if(Xind(i, k - 1) != Xind(i, k)) {
                    for(j = k; j < t.n_elem; j++) {
                        nXcum(j, Xind(i, k))++;
                    }
                }   
                nXsum(k, Xind(i, k))++;
            }
            
            // print outputs if required
            if(print == 1) {
                Rprintf("nXsum:\n");
                for(int r = 0; r < t.n_elem; r++) {
                    for(int c = 0; c < 12; c++) {
                        Rprintf("%d ", nXsum(r, c));
                    }
                    Rprintf("\n");
                }
                
                Rprintf("nXcum:\n");
                for(int r = 0; r < t.n_elem; r++) {
                    for(int c = 1; c < 12; c++) {
                        Rprintf("%d ", nXcum(r, c));
                    }
                    Rprintf("\n");
                }
            }
            
            // reorder old susceptibles
            if(Xind(i, t.n_elem - 1) == 0) {
                if(i < (Ninf - 1)) {
                    for(j = (i + 1); j < Ninf; j++) {
                        for(k = 0; k < t.n_elem; k++) {
                            Xind(j - 1, k) = Xind(j, k);
                        }
                    }
                }
                for(j = 0; j < t.n_elem; j++) {
                    Xind(Ninf - 1, j) = 0;
                }
                // update infection counters
                Ninf--;
            } else {
                // increment counter
                i++;
            }
        }
        
        /////////////////////////////////////////////////////////////////
        ///////                  SUSCEPTIBLES                     ///////
        /////////////////////////////////////////////////////////////////
        
//        Rprintf("sus\n");
        
        // update infection counters
        Nsus = Npop - Ninf_prev;
        
        // calculate filtering probabilities for remaining susceptibles
        i = Ninf;
        Ntemp = 0;
        while(Ntemp <= Nsus && i < Nposs) {
        
            // update transition probabilities
            updateTransProbs(&pXsum, &nXsum, &pars, Npop, pini);
            
            for(j = 0; j < t.n_elem; j++) {
                if(Xind(i, j) != 0) {
                    Rcpp::stop("mnahmnah\n");
                }
            }
            
            // remove individual i from states
            // and cumulative states to get correct
            // powers in subsequent updates
            nXsum(0, Xind(i, 0))--;
            for(k = 0; k < t.n_elem; k++) nXcum(k, Xind(i, 0))--;
            for(k = 1; k < t.n_elem; k++) {
                if(Xind(i, k - 1) != Xind(i, k)) {
                    for(j = k; j < t.n_elem; j++) nXcum(j, Xind(i, k))--;
                }
                nXsum(k, Xind(i, k))--;
            }

            // forwards filter
            forwardsFilter(i, pars(9), pars(10), &Xind, &pXind, &nXsum, &nXcum, &pXsum, &t, &D,  Npop, pini, &tempProb, print);
            
            // now extract probability of remaining susceptible
            psus = exp(pXind(i, t.n_elem - 1, 0));
            
            // take geometric random sample to give number of 
            // individuals who remain susceptible in this pass
            Ntemp += R::rgeom(1.0 - psus);
//            Rprintf("Ntemp = %d Nsus = %d\n", Ntemp, Nsus);
            
            // add the new infection
            Ntemp++;
            if(Ntemp < 0) Rcpp::stop("crap");
//            Rprintf("psus = %f Ntemp = %d Nsus = %d\n", psus, Ntemp, Nsus);
            
            // infect an individual if required
            if(Ntemp <= Nsus) {
                
                // update final time point filtering probs
                pXind(i, t.n_elem - 1, 0) = R_NegInf;
                
                // re-normalise probabilities
                norm = pXind(i, t.n_elem - 1, 0);
                for(k = 1; k < 12; k++) {
                    norm = (norm > pXind(i, t.n_elem - 1, k) ? norm:pXind(i, t.n_elem - 1, k));
                }
                norm1 = 0.0;
                for(k = 0; k < 12; k++) {
                    norm1 += exp(pXind(i, t.n_elem - 1, k) - norm);
                }
                norm += log(norm1);
                if(norm <= R_NegInf) {
                    Rcpp::stop("worrying backwards norm\n");
                }
                for(k = 0; k < 12; k++) {
                    pXind(i, t.n_elem - 1, k) -= norm;
                }
            
                // backwards sampler
                backwardsSampler (i, pars(9), pars(10), &pXind, &t, &Xind, &pXsum, &nXsum, &pXback, Npop, pini, print);
                
                // add individual i from states
                // and cumulative states to get correct
                // powers in subsequent updates
                nXsum(0, Xind(i, 0))++;
                for(k = 0; k < t.n_elem; k++) nXcum(k, Xind(i, 0))++;
                for(k = 1; k < t.n_elem; k++) {
                    if(Xind(i, k - 1) != Xind(i, k)) {
                        for(j = k; j < t.n_elem; j++) nXcum(j, Xind(i, k))++;
                    }   
                    nXsum(k, Xind(i, k))++;
                }
                
                // print outputs if required
                if(print == 1) {
                    Rprintf("nXsum:\n");
                    for(int r = 0; r < t.n_elem; r++) {
                        for(int c = 0; c < 12; c++) {
                            Rprintf("%d ", nXsum(r, c));
                        }
                        Rprintf("\n");
                    }
                    
                    Rprintf("nXcum:\n");
                    for(int r = 0; r < t.n_elem; r++) {
                        for(int c = 1; c < 12; c++) {
                            Rprintf("%d ", nXcum(r, c));
                        }
                        Rprintf("\n");
                    }
                }
                
                // update infective counts
                Ninf++;
                i++;
                
                // resize matrices if necessary
                if(Ninf == Nposs && Nposs < Npop) {
                    Rprintf("Resize\n");
                    Nposs = 2 * Nposs;
                    Nposs = (Nposs > Npop ? Npop:Nposs);
                    Xind.resize(Nposs, t.n_elem);
                    pXind.resize(Nposs, t.n_elem, 12);
                    for(k = Ninf; k < Nposs; k++) {
                        for(j = 0; j < t.n_elem; j++) {
                            Xind(k, j) = 0;
                        }
                    }
                }
            } else {
                // add individual i from states
                // and cumulative states to get correct
                // powers in subsequent updates
                nXsum(0, Xind(i, 0))++;
                for(k = 0; k < t.n_elem; k++) {
                    nXcum(k, Xind(i, 0))++;
                }
                for(k = 1; k < t.n_elem; k++) {
                    if(Xind(i, k - 1) != Xind(i, k)) {
                        for(j = k; j < t.n_elem; j++) {
                            nXcum(j, Xind(i, k))++;
                        }
                    }   
                    nXsum(k, Xind(i, k))++;
                }
            }
        }
//        Rprintf("Ninf = %d\n", Ninf);
        
        // update transition probabilities
        updateTransProbs(&pXsum, &nXsum, &pars, Npop, pini);
        
        // save output
        if(fixPars == 0) {
            for(int c = 0; c < npars; c++) {
                output(iter, c) = pars(c);
            }
        
            for(int r = 0; r < t.n_elem; r++) {
                for(int c = 0; c < 12; c++) {
                    output(iter, npars + r * 12 + c) = (outputCum == 0 || c == 0 ? nXsum(r, c):nXcum(r, c));
                }
            }
            output(iter, npars + t.n_elem * 12) = loglikelihood(pars(9), pars(10), &pXsum, &nXsum, &nXcum, &t, pini, Npop, &D);
        } else {        
            for(int r = 0; r < t.n_elem; r++) {
                for(int c = 0; c < 12; c++) {
                    output(iter, r * 12 + c) = (outputCum == 0 || c == 0 ? nXsum(r, c):nXcum(r, c));
                }
            }
            output(iter, t.n_elem * 12) = loglikelihood(pars(9), pars(10), &pXsum, &nXsum, &nXcum, &t, pini, Npop, &D);
        }   
                    
        // print some output to the screen for book-keeping
        if ((iter + 1) % 100 == 0) {	
            //calculate block run time
            timer.step("");
            Rcpp::NumericVector res(timer);
            
            if(fixPars == 0) {
                // update acceptance rate
                accrate = ((double) nacc) / ((double) 100);
                nacc = 0;
                
                Rprintf("i = %d acc = %.2f time = %.2f secs \n", iter + 1, accrate, (res[timer_cnt] / 1e9) - prev_time);
            } else {
                Rprintf("i = %d time = %.2f secs \n", iter + 1, (res[timer_cnt] / 1e9) - prev_time);
            }
            
            //reset timer and acceptance rate counter
            prev_time = res[timer_cnt] / 1e9;
            timer_cnt++;
        }
        
        // calculations for adaptive proposal
        if((iter + 1) >= 100 && fixPars == 0) {
            if((iter + 1) == 100) {
                // check for acceptances
                accrate = ((double) cumacc) / ((double) 100);
                if(accrate > 0.0) {
                	calcPost(iter, npars, &tempmn, &meanmat, &meanmat1, output, &propVar);
                } else {
                    arma::mat outlist (1, 1);
                    outlist(0, 0) = NA_REAL;
                    Rprintf("No initial acceptances.\n");
                    return(outlist);
                }
            } else {
            	adaptUpdate(iter, npars, &tempmn, &meanmat, &meanmat1, arma::conv_to<arma::vec>::from(output.row(iter)), &propVar);
            }
            propVarChol = cholArma(propVar * adaptscale, &cholScale);
        }
    }
    timer.step("");
    Rcpp::NumericVector res(timer);
    Rprintf("Final time = %.2f secs \n", res[timer_cnt] / 1e9);
            
    // return alphas
    return output;
}
