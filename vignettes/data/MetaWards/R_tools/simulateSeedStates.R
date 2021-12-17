## FFBS for early stage data
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(coda)
library(parallel)
library(lubridate)
library(stringr)

## extract command line arguments
args <- commandArgs(TRUE)
if(length(args) > 0) {
    stopifnot(length(args) == 2)
    startdate <- args[1]
    seeddate <- args[2]
} else {
    stop("No arguments")
}

## load seeds file
deaths <- readRDS("inputs/seedsInfo.rds")

# ## set dates
# startdate <- "2020-02-01"
# seeddate <- "2020-03-06"

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

## round any on boundaries to prevent overflow
pars <- mutate(pars, across(starts_with("p"), ~ifelse(. == 1, 0.999999, .))) %>%
    mutate(across(starts_with("p"), ~ifelse(. == 0, 0.000001, .)))

######################################################
#### get state distributions across multiple LADs ####
######################################################

## source iFFBS function
sourceCpp("R_tools/iFFBS.cpp")

## set seeding time at which to draw
seed_time <- as.numeric(ymd(seeddate) - ymd(startdate))

## set pini
pini <- 700 / (sum(deaths$popsize) * seed_time)

## extract LADS without deaths
nodeaths <- filter(deaths, deaths == 0)

## extract LADs with deaths
deaths <- filter(deaths, deaths > 0) %>%
    mutate(nlads = NA)

## group LADS by quartile according to population size
nodeaths <- mutate(nodeaths, quartile = ntile(popsize, 4))
nodeaths <- group_by(nodeaths, quartile) %>%
    summarise(med_popsize = round(median(popsize), 0), nlads = n()) %>%
    inner_join(nodeaths, by = "quartile") %>%
    arrange(lad) %>%
    group_by(quartile) %>%
    mutate(lad_within = 1:n()) %>%
    ungroup()

## bind some average LADs to death data
deaths <- select(nodeaths, -lad_within) %>%
    group_by(quartile) %>%
    slice(1) %>%
    ungroup() %>%
    select(LAD19CD, popsize = med_popsize, lad = quartile, data, deaths, nlads) %>%
    rbind(deaths)

## run simulations across multiple design points in parallel
seeds <- mclapply(1:nrow(pars), function(i, pars, deaths, seed_time, pini) {
    pars <- slice(pars, i)
    reps <- pars$repeats[1]
    pars <- select(pars, pH, pHD, pI1, pI1D, pI1H, pI2, pP, pE, pEP, nu, nuA) %>%
        unlist()
    ## generate seed distributions
    seeds <- pmap(list(deaths$data, deaths$popsize, deaths$lad, deaths$nlads), function(deaths, Npop, lad, nlads, pars, reps, pini) {
        ## set parameters
        D_prime <- deaths$cumDeaths
        D_prime <- c(D_prime[1], diff(D_prime))
        
        ## run iFFBS algorithm
        post1 <- iFFBS(D_prime, 1:length(D_prime), 
                       Npop, 2000, pini, pars = pars, fixPars = 1)
        
        ## set names and remove burnin
        colnames(post1) <- c(apply(expand.grid(c("S", "E", "A", "RA", "P", "I1", "I2", "RI", "DI", "H", "RH", "DH"), 1:length(D_prime)), 1, function(x) paste0(x[1], "_", as.numeric(x[2]))), "loglike")
        post1 <- as_tibble(post1) %>%
            select(ends_with(paste0("_", seed_time))) %>%
            slice(-(1:500))
        
        if(nrow(post1) < (ifelse(is.na(nlads), 1, nlads) * reps)) {
            stop("Not enough random samples to set seeds")
        }
        
        if(is.na(nlads)) {
            post <- slice_sample(post1, n = reps) %>%
                mutate(rep = 1:n(), lad = lad, lad_within = NA)
        } else {
            post <- slice_sample(post1, n = reps * nlads) %>%
                mutate(rep = rep(1:reps, each = nlads), lad = lad, 
                    lad_within = rep(1:nlads, times = reps))
        }
        post
    }, pars = pars, reps = reps, pini = pini)
    seeds <- bind_rows(seeds)
    colnames(seeds) <- gsub(paste0("_", seed_time), "", colnames(seeds))
    seeds
}, pars = pars, deaths = deaths, seed_time = seed_time, pini = pini, mc.cores = detectCores())
names(seeds) <- pars$output
seeds <- bind_rows(seeds, .id = "output")

## amend output names
seeds <- distinct(seeds, output, rep) %>%
    group_by(output) %>%
    count() %>%
    inner_join(seeds, by = "output") %>%
    mutate(output = ifelse(n > 1, paste0(output, "x", str_pad(rep, 3, pad = "0")), output)) %>%
    select(!c(n, rep)) %>%
    ungroup()

## join back to all LADs
stopifnot(
    nrow(
        select(nodeaths, quartile, lad, lad_within) %>%
        anti_join(seeds, by = c("quartile" = "lad", "lad_within" = "lad_within"))
    ) == 0
)
seeds <- left_join(
        seeds, 
        select(nodeaths, quartile, lad, lad_within), 
        by = c("lad" = "quartile", "lad_within" = "lad_within")
    ) %>%
    select(!lad_within) %>%
    mutate(lad = ifelse(is.na(lad.y), lad, lad.y)) %>%
    select(!lad.y) %>%
    arrange(output, lad)
stopifnot(
    group_by(seeds, output) %>%
    mutate(lad1 = 1:n()) %>%
    mutate(diff = lad - lad1) %>%
    ungroup() %>%
    pluck("diff") %>%
    {all(. == 0)}
)

## reorder to match ncov.json
seeds <- select(seeds, output, LAD = lad, E, P, I1, I2, RI, DI, A, RA, H, RH, DH)

## remove zeros
seeds <- seeds[rowSums(seeds[, -c(1, 2)]) > 0, ]

## check all outputs have at least one source of infection
## across LADS
stopifnot(
    all(
        rowwise(seeds) %>%
        mutate(ninf = sum(c_across(starts_with(c("E", "A", "P", "I1", "I2"))))) %>%
        ungroup() %>%
        group_by(output) %>%
        summarise(ninf = sum(ninf), .groups = "drop") %>%
        pluck("ninf") > 0))

## save seeding file
write_csv(seeds, "inputs/seeds.csv", col_names = FALSE)
