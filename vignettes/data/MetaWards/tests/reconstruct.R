## load libraries
library(Rcpp)

## write Rcpp function to reconstruct counts from incidence
cppFunction('IntegerMatrix reconstruct(IntegerVector Einc, IntegerVector Iinc, IntegerVector Rinc,
    IntegerVector Dinc, IntegerVector IAinc, IntegerVector RAinc, IntegerVector IHinc, IntegerVector RHinc, 
    IntegerVector DHinc, IntegerVector ICinc, IntegerVector RCinc, IntegerVector DCinc) {
    
    // extract sizes
    int n = Einc.size();
    
    // set up output matrix
    IntegerMatrix output(n, 17);
    
    // reconstruct counts
    int E = 0, I = 0, R = 0, D = 0, IA = 0, RA = 0, IH = 0;
    int RH = 0, DH = 0, IC = 0, RC = 0, DC = 0;
    for(int i = 0; i < n; i++) {
    
        E += Einc[i] - Iinc[i] - IAinc[i];
        output(i, 0) = Einc[i];
        output(i, 1) = E;
        
        I += Iinc[i] - IHinc[i] - Rinc[i] - Dinc[i];
        output(i, 2) = Iinc[i];
        output(i, 3) = I;
        
        R += Rinc[i];
        output(i, 4) = R;
        
        D += Dinc[i];
        output(i, 5) = D;
        
        IA += IAinc[i] - RAinc[i];
        output(i, 6) = IAinc[i];
        output(i, 7) = IA;
        
        RA += RAinc[i];
        output(i, 8) = RA;
        
        IH += IHinc[i] - ICinc[i] - RHinc[i] - DHinc[i];
        output(i, 9) = IHinc[i];
        output(i, 10) = IH;
        
        RH += RHinc[i];
        output(i, 11) = RH;
        
        DH += DHinc[i];
        output(i, 12) = DH;
        
        IC += ICinc[i] - RCinc[i] - DCinc[i];
        output(i, 13) = ICinc[i];
        output(i, 14) = IC;
        
        RC += RCinc[i];
        output(i, 15) = RC;
        
        DC += DCinc[i];
        output(i, 16) = DC;
    }
    
    // return counts
    return(output);
}')
