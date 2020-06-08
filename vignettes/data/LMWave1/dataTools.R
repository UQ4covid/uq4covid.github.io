## load libraries
library(tidyverse)
library(data.table)

## function to convert samples to correct format for MetaWards
convertToMetaWards <- function(pars, repeats) {
    disease <- tibble(`beta[2]` = pars$r_zero / pars$infectious_time)
    disease$`beta[3]` <- disease$`beta[2]`
    disease$`progress[1]` <- 1 / pars$incubation_time
    disease$`progress[2]` <- 1
    disease$`progress[3]` <- 1 / (pars$infectious_time - 1)
    disease$`.scale_rate[1]` <- pars$lock_1_restrict
    disease$`.scale_rate[2]` <- pars$lock_2_release
    disease$repeats <- repeats
    disease$par <- pars$par
    disease
}

## function to create fingerprint
createFingerprint <- function(disease) {
    
    ## check input
    stopifnot(is.data.frame(disease))
    
    ## extract repeats
    disease$fingerprint <- select(disease, -repeats, -par) %>%
        # mutate_all(round, digits = 6) %>%
        mutate_all(as.character) %>%
        mutate_all(function(x) {
            if(all(nchar(x) == 1)) {
                x <- paste0(x, ".0")
            }
            gsub("\\.", "i", x)
        }) %>%
        apply(1, paste0, collapse = "v")
    
    ## add replicates
    disease <- mutate(disease, repeats = map(repeats, ~1:.)) %>%
        unnest(cols = repeats) %>%
        mutate(fingerprintRep = paste0(fingerprint, "x", repeats)) %>%
        rename(`repeat` = repeats)
    
    ## return amended data set
    disease
}

## function to correct runs from MetaWards
readRuns <- function(main_dir = "raw_outputs", maxday) {
    
    ## check inputs
    checkInput(main_dir, c("vector", "character"), 1)
    checkInput(maxday, c("vector", "numeric"), 1, int = TRUE, gte = 0)
    
    ## read outputs
    dir <- paste0(main_dir, "/results.csv.bz2")
    sims <- fread(dir)
    
    sims <- sims %>%
        group_by(fingerprint, `repeat`) %>%
        nest() %>%
        mutate(data = map(data, function(x, maxday) {
            
            ## complete time points
            x <- complete(x, day = 0:maxday) %>%
                select(day, S, E, I, R) %>%
                mutate(D = 0)
            
            ## correct for initial conditions
            stopifnot((x$I[x$day == 1] + x$R[x$day == 1]) == 0)
            x$D[x$day == 1] <- x$E[x$day == 2]
            x$R[x$day == 1] <- x$D[x$day == 1]
            x$S[x$day == 1] <- x$S[x$day == 2]
            x <- mutate(x, D = lag(S, default = 0) - S)
            x$D[x$day == 0] <- 0
            
            ## update R and E classes
            x <- mutate(x, R = R - D) %>%
                mutate(E = E + D)
            
            ## extract INFECTION incidence 
            x <- rename(x, Einc = D) %>%
                mutate(Iinc = Einc - E + lag(E, default = 0)) %>%
                mutate(Rinc = R - lag(R, default = 0))
            
            ## generate cumulative counts
            x[is.na(x)] <- 0
            x <- mutate(x, Ecum = cumsum(Einc), Icum = cumsum(Iinc), Rcum = cumsum(Rinc))
            
            ## simple checks (necessary but not sufficient)
            stopifnot(all(x$Ecum >= x$Icum))
            stopifnot(all(x$Icum >= x$Rcum))
            stopifnot(all(x$Rcum == x$R))
            
            ## check counts
            stopifnot(all(
                select(x, S, E, I, R) %>%
                    as.matrix() %>%
                    apply(1, sum) == (x$S[1] + x$E[1])))
            x <- select(x, -Rcum) %>%
                filter(day > 0)
            x
        }, maxday = maxday)) %>%
        unnest(cols = data) %>%
        ungroup()
    sims
}
