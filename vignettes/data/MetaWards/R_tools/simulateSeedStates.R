## FFBS for early stage data
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(coda)
library(parallel)
library(lubridate)
library(stringr)

## load seeds file
deaths <- readRDS("inputs/seedsInfo.rds")

## set startdate (dmy)
startdate <- "06-03-2020"

## NGM for single population (order: E, A, P, I1, I2)
NGM <- function(R0 = NA, nu = NA, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2) {
  
    ## check that exactly one of R0 or nu is specified
    stopifnot(length(R0) == 1 & length(nu) == 1 & length(nuA) == 1)
    stopifnot(!is.na(R0) | !is.na(nu))
    stopifnot(is.na(R0) | is.na(nu))
    
    ## check inputs
    stopifnot(!missing(N) & !missing(S0))
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
    stopifnot(length(N) == 1 & length(S0) == 1)
    stopifnot(S0 <= N & S0 >= 1)
    
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
    F <- matrix(0, 5, 5)
    
    ## add A to E
    F[1, 2] <- S0 * nuA / N
    ## add P to E
    F[1, 3] <- S0 / N
    ## add I1 to E
    F[1, 4] <- S0 / N
    ## add I2 to E
    F[1, 5] <- S0 / N
    
    ## set up empty V matrix
    V <- matrix(0, 5, 5)
    
    ## add E to E terms
    V[1, 1] <- gammaE
    
    ## add E to A terms
    V[2, 1] <- -(1 - pEP) * gammaE
    ## add A to A terms
    V[2, 2] <- gammaA
    
    ## add E to P terms
    V[3, 1] <- -pEP * gammaE
    ## add P to P terms
    V[3, 3] <- gammaP
    
    ## add P to I1 terms
    V[4, 3] <- -gammaP
    ## add I1 to I1 terms
    V[4, 4] <- gammaI1
    
    ## add I1 to I2 terms
    V[5, 4] <- -(1 - pI1H - pI1D) * gammaI1
    ## add I2 to I2 terms
    V[5, 5] <- gammaI2
    
    ## remove pathways that don't exist
    rem <- NULL
    if(!pathA) rem <- 2
    if(!pathI1) rem <- c(rem, 3:5)
    if(!pathI2) rem <- c(rem, 5)
    if(!is.null(rem)) {
        rem <- unique(rem)
        F <- F[-rem, -rem]
        V <- V[-rem, -rem]
    }
    
    if(!is.na(R0)) {
        ## calculate NGM
        K <- F %*% solve(V)
        
        ## extract eigenvalues
        eig <- eigen(K)$values
        if(is.complex(eig)) {
            ## extract real eigenvalues
            eig <- as.numeric(eig[Im(eig) == 0])
            stopifnot(length(eig) >= 1)
        }
        
        ## calculate and return nu
        nu <- R0 / max(eig)
        return(list(nu = nu, K = nu * K))
    } else {        
        ## calculate NGM
        K <- nu * F %*% solve(V)
        
        ## extract eigenvalues
        eig <- eigen(K)$values
        if(is.complex(eig)) {
            ## extract real eigenvalues
            eig <- as.numeric(eig[Im(eig) == 0])
            stopifnot(length(eig) >= 1)
        }
        
        ## calculate and return R0
        R0 <- max(eig)
        return(list(R0 = R0, K = K))
    }
}

## load inputs
disease <- read_delim("inputs/disease.dat", delim = " ")
inputs <- readRDS("inputs/inputs.rds")
age_probs <- read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2

## read in parameters and calculate weighted average
## of transition parameters across age-classes
pars <- select(disease, ends_with(paste0("_", 1:8)), output) %>%
    pivot_longer(!output, names_to = "var", values_to = "value") %>%
    separate(var, c("var", "agegrp"), sep = "_") %>%
    mutate(agegrp = as.numeric(agegrp)) %>%
    mutate(ageprob = age_probs[agegrp]) %>%
    mutate(value = ageprob * value) %>%
    group_by(output, var) %>%
    summarise(value = sum(value), .groups = "drop") %>%
    pivot_wider(names_from = var, values_from = value)
colnames(pars) <- gsub("\\.", "", colnames(pars))

## rescale transmission parameter to reflect single population
pars <- pars %>%
    select(pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD, output) %>%
    inner_join(
        mutate(inputs, gammaE = 1 / TE, gammaA = 1 / (TP + TI1 + TI2), gammaP = 1 / TP, gammaI1 = 1 / TI1, gammaI2 = 1 / TI2), 
    by = "output") %>%
    mutate(nu = pmap_dbl(list(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2), 
        function(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2, S0, N) {
            NGM(R0 = R0, nu = NA, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2)$nu
        }, S0 = 149999, N = 150000)) %>%
    select(output, repeats, nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD) 

######################################################
#### get state distributions across multiple LADs ####
######################################################

## source iFFBS function
sourceCpp("R_tools/iFFBS.cpp")

## set seeding time at which to draw
seed_time <- as.numeric(dmy(startdate) - dmy("01-01-2020"))

## extract LADs with some deaths
deaths1 <- filter(deaths, deaths > 0) #%>%
    #rbind(filter(deaths, deaths == 0) %>% mutate(popsize = round(median(popsize))) %>% slice(1))

## run simulations across multiple design points in parallel
seeds <- mclapply(1:nrow(pars), function(i, pars, deaths, seed_time) {
    pars <- slice(pars, i)
    reps <- pars$repeats[1]
    pars <- select(pars, pH, pHD, pI1, pI1D, pI1H, pI2, pP, pE, pEP, nu, nuA) %>%
        unlist()
    ## generate seed distributions
    seeds <- map2(deaths$data, deaths$popsize, function(deaths, Npop, pars, reps) {
        ## set parameters
        D_prime <- deaths$cumDeaths
        D_prime <- c(D_prime[1], diff(D_prime))
        pini <- 20 / (Npop * length(D_prime))
        
        npost <- 0
        cycle <- 1
        post <- NULL
        if(npost < reps & cycle < 3) {
            ## run iFFBS algorithm
            post1 <- iFFBS(D_prime, 1:length(D_prime), 
                Npop, 2000, pini, pars = pars, fixPars = 1, outputCum = 1)
    
            colnames(post1) <- c(apply(expand.grid(c("S", "E", "A", "RA", "P", "I1", "I2", "RI", "DI", "H", "RH", "DH"), 1:length(D_prime)), 1, function(x) paste0(x[1], "_", as.numeric(x[2]))), "loglike")
            post1 <- as_tibble(post1) %>%
                select(ends_with(paste0("_", seed_time))) %>%
                slice(-(1:500))
            ## extract elements with at least one seed
            post1 <- rowwise(post1) %>%
                mutate(ninf = sum(c_across(starts_with(c("E_", "A_", "P_", "I1_", "I2_"))))) %>%
                ungroup() %>%
                filter(ninf > 0) %>%
                select(!ninf)
            post <- rbind(post, post1)
            npost <- nrow(post)
            cycle <- cycle + 1
        }   
        if(npost < reps & sum(D_prime) > 0) {
            stop(paste("Can't find enough non-zero states to sample from: npost = ", npost, ", nreps = ", reps))
        }
        post <- slice_sample(post, n = reps) %>%
            mutate(rep = 1:n())
        post
    }, pars = pars, reps = reps)
    names(seeds) <- deaths$lad
    seeds <- bind_rows(seeds, .id = "LAD")
    colnames(seeds) <- gsub(paste0("_", seed_time), "", colnames(seeds))
    seeds
}, pars = pars, deaths = deaths1, seed_time = seed_time, mc.cores = detectCores())
names(seeds) <- pars$output
seeds <- bind_rows(seeds, .id = "output")

## amend output names
seeds <- distinct(seeds, output, rep) %>%
    group_by(output) %>%
    count() %>%
    inner_join(seeds, by = "output") %>%
    mutate(output = ifelse(n > 1, paste0(output, "x", str_pad(rep, 3, pad = "0")), output)) %>%
    select(!c(n, rep))

## reorder to match ncov.json
seeds <- select(seeds, output, LAD, "E", "P", "I1", "I2", "RI", "DI", "A", "RA", "H", "RH", "DH")

## save seeding file
write_csv(seeds, "inputs/seeds.csv", col_names = FALSE)
