## load libraries
library(tidyverse)
library(Rcpp)

## create output directory
dir.create("outputs")

## NGM (order: E, A, P, I1, I2)
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

## read in parameters
pars <- read_delim("../../inputs/disease.dat", delim = " ") %>%
    select(ends_with("_8"), nu = `beta[1]`, nuA = `beta[6]` / `beta[1]`, output)
inputs <- readRDS("../../inputs/inputs.rds")
colnames(pars) <- gsub("_8", "", colnames(pars))
colnames(pars) <- gsub("\\.", "", colnames(pars))
pars <- pars %>%
    select(nuold = nu, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD, output) %>%
    inner_join(
        mutate(inputs, gammaE = 1 / TE, gammaA = 1 / (TP + TI1 + TI2), gammaP = 1 / TP, gammaI1 = 1 / TI1, gammaI2 = 1 / TI2), 
    by = "output") %>%
    mutate(nu = pmap_dbl(list(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2), 
        function(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2, S0, N) {
            NGM(R0 = R0, nu = NA, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2)$nu
        }, S0 = 149, N = 150))

## plot old scaling against new
p <- ggplot(pars) +
    geom_point(aes(x = nuold, y = nu)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")
ggsave("outputs/nuscaling.pdf", p)

## extract parameters for simulation   
pars <- select(slice(pars, 6), nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD) %>%
    unlist()
pars["pI2"] <- 0.5
pars["pHD"] <- 0.05

## set initial counts
u <- numeric(12)
u[1] <- 149
u[2] <- 1

## set seed
set.seed(4578)

## try discrete-time model
sourceCpp("discreteStochModel.cpp")
disSims <- list()
for(i in 1:50) {
    disSims[[i]] <- discreteStochModel(pars, 150, u)
}
stageNms <- c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
disSims <- map(disSims, ~as_tibble(.)) %>%
    bind_rows(.id = "rep") %>%
    set_names(c("rep", "t", stageNms)) %>%
    mutate(t = t + 1)

## extract simulation closest to median
medRep <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        median = median(n),
        .groups = "drop"
    ) %>%
    inner_join(
        pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n"),
        by = c("t", "var")
    ) %>%
    mutate(diff = abs(median - n)) %>%
    group_by(rep) %>%
    summarise(diff = sum(diff), .groups = "drop") %>%
    arrange(diff) %>%
    slice(2) %>%
    inner_join(disSims, by = "rep") %>%
    select(!c(diff, rep))
    
## plot replicates
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = median)) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(medRep, !t, names_to = "var", values_to = "n"),
            col = "red", linetype = "dashed"
        ) +
        facet_wrap(~var) +
        xlab("Days") + 
        ylab("Counts")
ggsave("outputs/sims.pdf", p)

saveRDS(medRep, "outputs/disSims.rds")
saveRDS(pars, "outputs/pars.rds")

