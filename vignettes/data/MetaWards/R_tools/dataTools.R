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

## @knitr maximin
## function to generate maximin samples given an arbitrary FMM
FMMmaximin <- function(model, nsamp, nseed = 100000) {
    
    ## check inputs and dependencies
    require(mclust)
    require(fields)
  
    stopifnot(!missing(model) & !missing(nsamp))
    stopifnot(class(model)[1] == "densityMclust")
    stopifnot(is.numeric(nsamp) & round(nsamp) == nsamp & nsamp > 1)
    stopifnot(is.numeric(nseed) & round(nseed) == nseed & nseed > 1 & nsamp < nseed)
    
    ## produce large number of samples from model
    sims <- sim(model$modelName, model$parameters, nseed)[, -1]
    
    ## rescale
    simsnorm <- apply(sims, 2, function(x) (x - min(x)) / diff(range(x)))
    
    ## sample initial choice at random
    simind <- sample.int(nrow(sims), 1)
    simrem <- c(1:nrow(sims))[-simind]
    
    ## choose other points using maximin
    while(length(simind) < nsamp) {
        dists <- rdist(simsnorm[simind, , drop = FALSE], simsnorm[-simind, , drop = FALSE])
        dists <- apply(dists, 2, min)
        dists <- which(dists == max(dists))
        simind <- c(simind, simrem[dists[1]])
        simrem <- simrem[-dists[1]]
    }
    
    ## return design
    sims[simind, ]
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

## @knitr NGM
## NGM (order: E, A, P, I1, I2)
NGM <- function(R0 = NA, nu = NA, C, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2) {
  
    ## check that exactly one of R0 or nu is specified
    stopifnot(!is.na(R0) | !is.na(nu))
    stopifnot(is.na(R0) | is.na(nu))
    
    ## check inputs
    stopifnot(!missing(C) & !missing(N) & !missing(S0))
    if(missing(nuA)) nuA <- 0
    if(missing(gammaE)) gammaE <- 0
    if(missing(pEP)) pEP <- 0
    if(missing(gammaA)) gammaA <- 0
    if(missing(gammaP)) gammaP <- 0
    if(missing(gammaI1)) gammaI1 <- 0
    if(missing(pI1H)) pI1H <- 0
    if(missing(pI1D)) pI1D <- 0
    if(missing(gammaI2)) gammaI2 <- 0
    stopifnot(all(c(nuA, gammaE, pEP, gammaA, 
                    gammaP, gammaI1, pI1H, pI1D, gammaI2) >= 0))
    
    ## extract lengths
    nage <- length(N)
    stopifnot(all(dim(C) == nage))
    stopifnot(length(S0) == nage & all(S0 <= N))
    
    ## check for empty pathways
    stopifnot(gammaE > 0)
    pathA <- ifelse(pEP == 1, FALSE, TRUE)
    if(pathA) stopifnot(gammaA > 0)
    pathI1 <- ifelse(pEP == 0, FALSE, TRUE)
    if(pathI1) stopifnot(gammaP > 0 & gammaI1 > 0)
    pathI2 <- ifelse((1 - pI1H - pI1D) == 0 | !pathI1, FALSE, TRUE)
    if(pathI2) stopifnot(gammaI2 > 0)
    
    ## set up empty F matrix
    ## (order: E, A, P, I1, I2)
    F <- matrix(0, nage * 5, nage * 5)
    
    ## add A to E
    F[1:nage, (nage + 1):(2 * nage)] <- t(t(nuA * S0 * C) / ifelse(N == 0, 1, N))
    ## add P to E
    F[1:nage, (2 * nage + 1):(3 * nage)] <- t(t(S0 * C) / ifelse(N == 0, 1, N))
    ## add I1 to E
    F[1:nage, (3 * nage + 1):(4 * nage)] <- t(t(S0 * C) / ifelse(N == 0, 1, N))
    ## add I2 to E
    F[1:nage, (4 * nage + 1):(5 * nage)] <- t(t(S0 * C) / ifelse(N == 0, 1, N))
    
    ## set up empty V matrix
    V <- matrix(0, nage * 5, nage * 5)
    
    ## add E to E terms
    diag(V)[1:nage] <- gammaE
    
    ## add E to A terms
    V[(nage + 1):(2 * nage), 1:nage] <- diag(nage) * -(1 - pEP) * gammaE
    ## add A to A terms
    V[(nage + 1):(2 * nage), (nage + 1):(2 * nage)] <- diag(nage) * gammaA
    
    ## add E to P terms
    V[(2 * nage + 1):(3 * nage), 1:nage] <- diag(nage) * -pEP * gammaE
    ## add P to P terms
    V[(2 * nage + 1):(3 * nage), (2 * nage + 1):(3 * nage)] <- diag(nage) * gammaP
    
    ## add P to I1 terms
    V[(3 * nage + 1):(4 * nage), (2 * nage + 1):(3 * nage)] <- diag(nage) * -gammaP
    ## add I1 to I1 terms
    V[(3 * nage + 1):(4 * nage), (3 * nage + 1):(4 * nage)] <- diag(nage) * gammaI1
    
    ## add I1 to I2 terms
    V[(4 * nage + 1):(5 * nage), (3 * nage + 1):(4 * nage)] <- diag(nage) * -(1 - pI1H - pI1D) * gammaI1
    ## add I2 to I2 terms
    V[(4 * nage + 1):(5 * nage), (4 * nage + 1):(5 * nage)] <- diag(nage) * gammaI2
    
    ## remove pathways that don't exist
    rem <- NULL
    if(!pathA) rem <- (nage + 1):(2 * nage)
    if(!pathI1) rem <- c(rem, (2 * nage + 1):(5 * nage))
    if(!pathI2) rem <- c(rem, (4 * nage + 1):(5 * nage))
    if(!is.null(rem)) {
        rem <- unique(rem)
        F <- F[-rem, -rem]
        V <- V[-rem, -rem]
    }
    
    if(!is.na(R0)) {
        ## calculate NGM
        K <- F %*% solve(V)
        
        ## calculate and return nu
        nu <- R0 / max(eigen(K)$values)
        return(list(nu = nu, K = nu * K))
    } else {        
        ## calculate NGM
        K <- nu * F %*% solve(V)
        
        ## calculate and return nu
        R0 <- max(eigen(K)$values)
        return(list(R0 = R0, K = K))
    }
}
