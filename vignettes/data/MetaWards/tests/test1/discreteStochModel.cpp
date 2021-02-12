#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
IntegerMatrix discreteStochModel(NumericVector pars, int tstart, int tstop, IntegerVector u) {
    
    int i = 0;
  
    // extract parameters
    double beta = pars[0];
    double gammaE = pars[1];
    double gammaP = pars[2];
    double gammaI1 = pars[3];
      
    // set up output
    IntegerVector u1(u.size());
    for(i = 0; i < u.size(); i++) u1[i] = u[i];
    
    IntegerMatrix out (tstop - tstart + 1, u.size() + 1);
    out(0, 0) = 0;
    for(i = 0; i < u.size(); i++) out(0, i + 1) = u[i];
    
    NumericVector rates(4);
    double N = sum(u);
    rates[0] = beta * (u[2] + u[3]) / N;
    rates[1] = gammaE;
    rates[2] = gammaP;
    rates[3] = gammaI1;
    double totrate = sum(rates);
    double prob = 0.0;
    tstart++;
    
    while(tstart < tstop) {

        // discrete-time model

        if(totrate > 0.0) {
            // SE
            prob = 1 - exp(-rates[0]);
            i = R::rbinom(u[0], prob);
            u1[0] -= i;
            u1[1] += i;
            
            // EP
            prob = 1 - exp(-rates[1]);
            i = R::rbinom(u[1], prob);
            u1[1] -= i;
            u1[2] += i;
            
            // PI1
            prob = 1 - exp(-rates[2]);
            i = R::rbinom(u[2], prob);
            u1[2] -= i;
            u1[3] += i;
            
            // I1DI
            prob = 1 - exp(-rates[3]);
            i = R::rbinom(u[3], prob);
            u1[3] -= i;
            u1[4] += i;
        }
        
        // record output
        out(tstart, 0) = tstart;
        for(i = 0; i < u.size(); i++) {
            out(tstart, i + 1) = u1[i];
            u[i] = u1[i];
        }
        tstart++;
        
        // update rates
        rates[0] = beta * (u[2] + u[3]) / N;
        rates[1] = gammaE;
        rates[2] = gammaP;
        rates[3] = gammaI1;
        totrate = sum(rates);
    }
    
    return out;
}
