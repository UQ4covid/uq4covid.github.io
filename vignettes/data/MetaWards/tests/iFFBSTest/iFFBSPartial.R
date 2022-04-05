## FFBS for early stage data
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(coda)
library(furrr)
library(lubridate)

## set path to MetaWardsData
path <- paste0("../../../../../../MetaWardsData/model_data/2011to2019Data/")

## load workers and players to obtain population sizes
workers <- read_delim(paste0(path, "WorkSize19.dat"), col_names = FALSE, delim = " ")
players <- read_delim(paste0(path, "PlaySize19.dat"), col_names = FALSE, delim = " ")
popsize <- full_join(workers, players, by = "X1") %>%
    replace_na(list(X2.x = 0, X2.y = 0)) %>%
    mutate(popsize = X2.x + X2.y) %>%
    select(-X2.x, -X2.y)

## now read in ward to LAD lookup
Ward19Lookup <- read_csv(paste0(path, "Ward19_Lookup.csv"))

## extract LADS
LADs <- unique(Ward19Lookup$LAD19CD)

## join to get population size by LAD
popsize <- inner_join(popsize, Ward19Lookup, by = c("X1" = "FID")) %>%
    group_by(LAD19CD) %>%
    summarise(popsize = sum(popsize), .groups = "drop")

## read in data, originally from here: 
## https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=cumDeathsByDeathDate&format=csv
deaths <- read_csv("../../data/seedDeaths/ltla_2021-05-18.csv", guess_max = 30000) %>%
    filter(date <= "2020-03-13") %>%
    select(areaCode, date, cumDeathsByDeathDate) %>%
    arrange(areaCode, date) %>%
    group_by(areaCode) %>%
    mutate(deaths = cumDeathsByDeathDate - lag(cumDeathsByDeathDate)) %>%
    mutate(deaths = ifelse(row_number() == 1, cumDeathsByDeathDate, deaths)) %>%
    ungroup() %>%
    complete(areaCode, date = seq(as.Date("2020-01-01"), as.Date("2020-03-13"), by = 1),
             fill = list(deaths = 0)) %>%
    select(-cumDeathsByDeathDate) %>%
    group_by(areaCode) %>%
    nest() %>%
    inner_join(popsize, by = c("areaCode" = "LAD19CD")) %>%
    mutate(totdeaths = map_dbl(data, ~sum(.$deaths)))

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
disease <- read_delim("../../inputs/disease.dat", delim = " ")
inputs <- readRDS("../../inputs/inputs.rds")
age_probs <- read_csv("../../inputs/age_seeds.csv", col_names = FALSE)$X2

## read in parameters and check scaling
output <- list()
for(i in c(1, NA, 8)) {
    if(!is.na(i)) {
        pars <- select(disease, ends_with(paste0("_", i)), nu = `beta[1]`, output)
        colnames(pars) <- gsub(paste0("_", i), "", colnames(pars))
        colnames(pars) <- gsub("\\.", "", colnames(pars))
    } else {
        ## weighted average
        pars <- select(disease, ends_with(paste0("_", 1:8)), nu = `beta[1]`, output) %>%
            pivot_longer(!c(nu, output), names_to = "var", values_to = "value") %>%
            separate(var, c("var", "agegrp"), sep = "_") %>%
            mutate(agegrp = as.numeric(agegrp)) %>%
            mutate(ageprob = age_probs[agegrp]) %>%
            mutate(value = ageprob * value) %>%
            group_by(nu, output, var) %>%
            summarise(value = sum(value), .groups = "drop") %>%
            pivot_wider(names_from = var, values_from = value)
        colnames(pars) <- gsub("\\.", "", colnames(pars))
    }
    pars <- pars %>%
        select(nuold = nu, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD, output) %>%
        inner_join(
            mutate(inputs, gammaE = 1 / TE, gammaA = 1 / (TP + TI1 + TI2), gammaP = 1 / TP, gammaI1 = 1 / TI1, gammaI2 = 1 / TI2), 
        by = "output") %>%
        mutate(nu = pmap_dbl(list(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2), 
            function(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2, S0, N) {
                NGM(R0 = R0, nu = NA, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2)$nu
            }, S0 = 149999, N = 150000))
    output[[length(output) + 1]] <- pars
}
names(output) <- c("young", "weighted", "old")
output <- bind_rows(output, .id = "group")

## plot old scaling against new
p <- ggplot(output) +
    geom_point(aes(x = nuold, y = nu)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    facet_wrap(~group)
ggsave("outputs/nuscalingPartial.pdf", p)

#######################################################
#### run multi parameter test across multiple LADs ####
#######################################################

## set parallel cores
plan(multisession, workers = 8)
furrr_options(seed = 666)

## run simulations across multiple design points
## in parallel
seeds <- filter(output, output %in% output[1:9]) %>%
    select(output, group, pH, pHD, pI1, pI1D, pI1H, pI2, pP, pE, pEP, nu, nuA) %>%
    group_by(output, group) %>%
    nest() %>%
    mutate(seeds = map(data, ~{
        apply(., 1, function(pars, deaths) {
            ## generate seeds
            seeds <- future_map2(deaths$data, deaths$popsize, function(deaths, Npop, pars) {
                ## source FF function
                sourceCpp("iFFBSPartialMCMC.cpp")
                
                ## set parameters
                D_prime <- deaths$deaths
                pini <- 1 / (Npop * length(D_prime))
                
                ## run iFFBS algorithm
                post <- iFFBS(D_prime, 1:length(D_prime), 
                    Npop, 2000, pini, pars = pars, fixPars = 1, outputCum = 1)

                colnames(post) <- c(apply(expand.grid(c("S", "E", "A", "RA", "P", "I1", "I2", "RI", "DI", "H", "RH", "DH"), 1:length(D_prime)), 1, function(x) paste0(x[1], "_", as.numeric(x[2]))), "loglike")
                post
            }, pars = pars)
            seeds
        }, deaths = deaths[1:10, ])
    }))

## extract states at final time point
temp <- map(seeds$seeds, ~{
        map(.[[1]], ~{
            as_tibble(.) %>%
            select(ends_with("_73")) %>%
            slice(-(1:500))
        }) %>%
        set_names(deaths$areaCode[1:10]) %>%
        bind_rows(.id = "LAD")
    }) %>%
    set_names(paste0(seeds$output, "_", seeds$group)) %>%
    bind_rows(.id = "output") %>%
    separate(output, c("output", "group"), sep = "_")

## plot cumE at time 73 across LADS
p <- ggplot(temp) +
    geom_violin(aes(x = group, y = E_73)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(LAD ~ output) +
    ylab("Cumulative E at day 73") + xlab("Group")
ggsave("outputs/cumE73.pdf", p, width = 15, height = 10)
    
## extract time of first introduction
firstTime <- map(seeds$seeds, ~{
        map(.[[1]], ~{
            as_tibble(.) %>%
                select(starts_with("E_")) %>%
                apply(1, function(x) {
                    which(x > 0)[1]
                }) %>%
                {tibble(day = .)} %>%
                mutate(iter = 1:n())
        }) %>%
        set_names(deaths$areaCode[1:10]) %>%
        bind_rows(.id = "LAD")
    }) %>%
    set_names(paste0(seeds$output, "_", seeds$group)) %>%
    bind_rows(.id = "output") %>%
    separate(output, c("output", "group"), sep = "_") %>%
    mutate(Date = dmy("01/01/2020") + day - 1)

## plot initial introduction
p <- ggplot(firstTime) +
    geom_violin(aes(x = Date, y = group)) +
    facet_grid(LAD ~ output) + ylab("Group") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("outputs/initialIntro.pdf", p, width = 15, height = 10)
        
## plot traces
pdf("outputs/traces.pdf")
group_by(firstTime, output, group) %>%
    nest() %>%
    mutate(data = map(data, ~{
        select(., iter, LAD, day) %>%
        pivot_wider(names_from = LAD, values_from = day) %>%
        select(-iter) %>%
        as.matrix() %>%
        as.mcmc()
    })) %>%
    {pmap(list(.$output, .$group, .$data), function(output, group, samples) {
        colnames(samples) <- paste0(colnames(samples), "_", output, "_", group)
        plot(samples, density = FALSE)
    })}
dev.off()

