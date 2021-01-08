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
      filter(alphaHC < -etaHC * ageM) %>%
      filter(alphaHR < -etaHR * ageM) %>%
      filter(alphaCR < -etaCR * ageM)
    
    ## return inputs
    input
}

## @knitr convertInputToDisease
## function to convert from input to disease space
convertInputToDisease <- function(input, C, ages) {
   
    require(dplyr)
    
    ## transmission parameter
    disease <- tibble(`.nu` = (input$r_zero / (max(eigen(C)$values) * input$infectious_time)))
    stopifnot(all(disease$`.nu` > 0 & disease$`.nu` < 1))
    
    ## progressions out of the E class
    disease$`.pE` <- 1 - exp(-1 / input$incubation_time)
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
    disease$`.pA` <- 1 - exp(-1 / input$infectious_time)
    
    ## progressions out of the I class
    disease$`.pI` <- 1 - exp(-1 / input$infectious_time)
    disease <- mutate(disease, 
        temp = map2(input$alphaIH, input$etaIH, function(alpha, eta, ages) {
            out <- exp(alpha + eta * ages) %>%
                as.matrix(nrow = 1) %>%
                as_tibble()
            colnames(out) <- paste0(".pIH_", 1:length(ages))
            out
        }, ages = ages)) %>%
        unnest(cols = temp)
    disease <- mutate(disease, 
        temp = map2(input$alphaIR, input$etaIR, function(alpha, eta, ages) {
            out <- exp(alpha + eta * ages) %>%
                as.matrix(nrow = 1) %>%
                as_tibble()
            colnames(out) <- paste0(".pIR_", 1:length(ages))
            out
        }, ages = ages)) %>%
        unnest(cols = temp)
    
    ## progressions out of the H class
    disease$`.pH` <- 1 - exp(-1 / input$hospital_time)
    disease <- mutate(disease, 
        temp = map2(input$alphaHC, input$etaHC, function(alpha, eta, ages) {
            out <- exp(alpha + eta * ages) %>%
                as.matrix(nrow = 1) %>%
                as_tibble()
            colnames(out) <- paste0(".pHC_", 1:length(ages))
            out
        }, ages = ages)) %>%
        unnest(cols = temp)
    disease <- mutate(disease, 
        temp = map2(input$alphaHR, input$etaHR, function(alpha, eta, ages) {
            out <- exp(alpha + eta * ages) %>%
                as.matrix(nrow = 1) %>%
                as_tibble()
            colnames(out) <- paste0(".pHR_", 1:length(ages))
            out
        }, ages = ages)) %>%
        unnest(cols = temp)
    
    ## progressions out of the C class
    disease$`.pC` <- 1 - exp(-1 / input$critical_time)
    disease <- mutate(disease, 
          temp = map2(input$alphaCR, input$etaCR, function(alpha, eta, ages) {
          out <- exp(alpha + eta * ages) %>%
              as.matrix(nrow = 1) %>%
              as_tibble()
          colnames(out) <- paste0(".pCR_", 1:length(ages))
          out
      }, ages = ages)) %>%
      unnest(cols = temp)
    
    ## lockdown scalings
    disease$`.lock_1_restrict` <- input$lock_1_restrict
    disease$`.lock_2_release` <- input$lock_2_release
    
    ## set up scaling for mixing matrix
    disease$`.GP_A` <- input$GP_A
    disease$`.GP_H` <- input$GP_H
    disease$`.GP_C` <- input$GP_C
    
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

