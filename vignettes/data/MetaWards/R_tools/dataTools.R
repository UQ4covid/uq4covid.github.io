## @knitr convertDesignToInput
## function to convert from design to input space
convertDesignToInput <- function(design, parRanges, scale = c("zero_one", "negone_one")) {
  
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
    
    ## return inputs
    input
}

## @knitr convertInputToDisease
## function to convert from input to disease space
convertInputToDisease <- function(input, C) {
   
    require(dplyr)
    
    ## transmission parameter
    disease <- tibble(`.nu` = (input$r_zero / (max(eigen(C)$values) * input$infectious_time)))
    
    ## progressions out of the E class
    ## (adjusting in the way described in the vignette for
    ## order of movers)
    disease$`.pE` <- 1 - exp(-1 / input$incubation_time)
    disease$`.pEA` <- input$pEA
    
    ## progressions out of the A class
    disease$`.pA` <- 1 - exp(-1 / input$infectious_time)
    
    ## progressions out of the I class
    ## (adjusting in the way described in the vignette for
    ## order of movers)
    disease$`.pI` <- 1 - exp(-1 / input$infectious_time)
    disease$`.pIH` <- input$pIH
    disease$`.pIR` <- (1 - input$pIH) * input$pIRprime
    
    ## progressions out of the H class
    ## (adjusting in the way described in the vignette for
    ## order of movers)
    disease$`.pH` <- 1 - exp(-1 / input$hospital_time)
    disease$`.pHC` <- input$pHC
    disease$`.pHR` <- (1 - input$pHC) * input$pHRprime
    
    ## progressions out of the C class
    ## (adjusting in the way described in the vignette for
    ## order of movers)
    disease$`.pC` <- 1 - exp(-1 / input$critical_time)
    disease$`.pCR` <- input$pCR
    
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
