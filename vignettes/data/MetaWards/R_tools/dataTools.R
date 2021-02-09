## @knitr convertDesignToInput
## function to convert from design to input space
convertDesignToInput <- function(design, parRanges, ageM, scale = c("zero_one", "negone_one")) {
  
    require(dplyr)
    require(tidyr)
  
    ## convert from design space to input space
    input <- mutate(design, ind = 1:n()) %>%
        gather(parameter, value, -output, -ind, -repeats) %>%
        left_join(parRanges, by = "parameter")
    if(scale[1] == "negone_one") {
        input <- mutate(input, value = value / 2 + 1)
    }
    input <- mutate(input, value = value * (upper - lower) + lower) %>%
        dplyr::select(ind, output, repeats, parameter, value) %>%
        spread(parameter, value) %>%
        arrange(ind) %>%
        dplyr::select(-ind)
    
    ## remove design points that don't adhere to necessary criteria
    input <- filter(input, alphaEA < -etaEA * ageM) %>%
      filter(alphaIH < -etaIH * ageM) %>%
      filter(alphaIR < -etaIR * ageM) %>%
      filter(alphaHR < -etaHR * ageM)
    
    ## return inputs
    input
}

## @knitr convertInputToDisease
## function to convert from input to disease space
convertInputToDisease <- function(input, C, ages) {
   
    require(dplyr)
    
    ## transmission parameter
    disease <- tibble(`.nu` = (input$r_zero / (max(eigen(C)$values) * input$infectious1_time)))
    stopifnot(all(disease$`.nu` > 0 & disease$`.nu` < 1))
    
    ## progressions out of the E class
    for(j in 1:length(ages)) {
        disease <- mutate(disease, ".pE_{j}" := 1 - exp(-1 / input$latent_time))
    }
    disease <- mutate(disease, 
        temp = map2(input$alphaEA, input$etaEA, function(alpha, eta, ages) {
            out <- exp(alpha + eta * ages) %>%
                as.matrix(nrow = 1) %>%
                as_tibble()
            colnames(out) <- paste0(".pEA_", 1:length(ages))
            out
        }, ages = ages)) %>%
        unnest(cols = temp)
    
    ## progressions out of the A class
    for(j in 1:length(ages)) {
        disease <- mutate(disease, ".pA_{j}" := 1 - exp(-1 / input$infectious1_time))
    }
    
    ## progressions out of the I1 class
    for(j in 1:length(ages)) {
        disease <- mutate(disease, ".pI1_{j}" := 1 - exp(-1 / input$infectious1_time))
    }
    disease <- mutate(disease, 
        temp = map2(input$alphaI1H, input$etaI1H, function(alpha, eta, ages) {
            out <- exp(alpha + eta * ages) %>%
                as.matrix(nrow = 1) %>%
                as_tibble()
            colnames(out) <- paste0(".pI1H_", 1:length(ages))
            out
        }, ages = ages)) %>%
        unnest(cols = temp)
    disease <- mutate(disease, 
        temp = map2(input$alphaI1I2, input$etaI1I2, function(alpha, eta, ages) {
            out <- exp(alpha + eta * ages) %>%
                as.matrix(nrow = 1) %>%
                as_tibble()
            colnames(out) <- paste0(".pI1I2_", 1:length(ages))
            out
        }, ages = ages)) %>%
        unnest(cols = temp)
    
    ## progressions out of the I2 class
    for(j in 1:length(ages)) {
        disease <- mutate(disease, ".pI2_{j}" := 1 - exp(-1 / input$infectious2_time))
    }
    
    ## progressions out of the H class
    for(j in 1:length(ages)) {
        disease <- mutate(disease, ".pH_{j}" := 1 - exp(-1 / input$hospital_time))
    }
    disease <- mutate(disease, 
        temp = map2(input$alphaHR, input$etaHR, function(alpha, eta, ages) {
            out <- exp(alpha + eta * ages) %>%
                as.matrix(nrow = 1) %>%
                as_tibble()
            colnames(out) <- paste0(".pHR_", 1:length(ages))
            out
        }, ages = ages)) %>%
        unnest(cols = temp)
    
    ## lockdown scalings
    disease$`.lock_1_restrict` <- input$lock_1_restrict
    disease$`.lock_2_release` <- input$lock_2_release
    
    ## set up scaling for mixing matrix
    disease$`.GP_A` <- input$GP_A
    disease$`.GP_H` <- input$GP_H
    
    ## finalise number of repeats
    disease$repeats <- input$repeats
    
    ## copy hash for each input
    disease$output <- input$output
    
    ## return disease file
    disease
}

## @knitr ensembleIDGen
## function for generating ensemble hash IDs (from Danny)
ensembleIDGen <- function(ensembleID = "a0", ensembleSize) {
    HexString <- c(as.character(0:9), letters)
    counter1 <- 1
    counter2 <- 1
    counter3 <- 1
    TotalGenerated <- 0
    tIDs <- rep(ensembleID, ensembleSize)
    while(counter1 <= 36 & TotalGenerated < ensembleSize){
      while(counter2 <= 36 & TotalGenerated < ensembleSize){
        while(counter3 <= 36 & TotalGenerated < ensembleSize){
          TotalGenerated <- TotalGenerated + 1
          tIDs[TotalGenerated] <- paste0(tIDs[TotalGenerated], HexString[counter1], HexString[counter2], HexString[counter3])
          counter3 <- counter3 + 1
        }
        counter2 <- counter2 + 1
        counter3 <- 1
      }
      counter1 <- counter1 + 1
      counter2 <- 1
      counter3 <- 1
    }
    stopifnot(!any(duplicated(tIDs)))
    tIDs
}

## @knitr reconstruct
## write Rcpp function to reconstruct counts from incidence
library(Rcpp)
cppFunction('IntegerMatrix reconstruct(IntegerVector Einc, IntegerVector Pinc, 
    IntegerVector I1inc, IntegerVector I2inc, IntegerVector RIinc, IntegerVector DIinc, 
    IntegerVector Ainc, IntegerVector RAinc, 
    IntegerVector Hinc, IntegerVector RHinc, IntegerVector DHinc) {
    
    // extract sizes
    int n = Einc.size();
    
    // set up output matrix
    IntegerMatrix output(n, 17);
    
    // reconstruct counts
    int E = 0, P = 0, I1 = 0, I2 = 0, RI = 0, DI = 0, A = 0, RA = 0, H = 0, RH = 0, DH = 0;
    for(int i = 0; i < n; i++) {
    
        E += Einc[i] - Pinc[i] - Ainc[i];
        output(i, 0) = Einc[i];
        output(i, 1) = E;
    
        P += Pinc[i] - I1inc[i];
        output(i, 2) = Pinc[i];
        output(i, 3) = P;
        
        I1 += I1inc[i] - I2inc[i] - Hinc[i] - DIinc[i];
        output(i, 4) = I1inc[i];
        output(i, 5) = I1;
        
        I2 += I2inc[i] - RIinc[i];
        output(i, 6) = I2inc[i];
        output(i, 7) = I2;
        
        RI += RIinc[i];
        output(i, 8) = RI;
        
        DI += DIinc[i];
        output(i, 9) = DI;
        
        A += Ainc[i] - RAinc[i];
        output(i, 10) = Ainc[i];
        output(i, 11) = A;
        
        RA += RAinc[i];
        output(i, 12) = RA;
        
        H += Hinc[i] - RHinc[i] - DHinc[i];
        output(i, 13) = Hinc[i];
        output(i, 14) = H;
        
        RH += RHinc[i];
        output(i, 15) = RH;
        
        DH += DHinc[i];
        output(i, 16) = DH;
    }
    
    // return counts
    return(output);
}')

