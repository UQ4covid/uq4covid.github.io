## load libraries
library(mclust)
library(MASS)
library(tidyverse)
library(coda)
library(GGally)
library(svglite)
library(parallel)

## read in data
IFR <- readRDS("IFR.rds")
hosp <- readRDS("hosp.rds")

## implements log-sum-exp trick
log_sum_exp <- function(x) {
    
    ## extract length
    n <- length(x)
    
    ## extract non-finite values
    ## (it's OK to remove these here because we
    ## assume they correspond to a likelihood of zero
    ## and are thus felt through 1 / n above)
    x <- x[is.finite(x)]
    if(length(x) == 0) {
      return(-Inf)
    }
    
    ## extract maximum of logged values
    mX <- max(x)
    
    ## return answer
    out <- mX + log(sum(exp(x - mX)))
    out <- out - log(n)
    out
}

## M-H algorithm
MH <- function(N, IFR, hosp, nind, nMC, LB = -20, UB = 0, nupdate = 100) {
    
    ## extract ageMid
    ageMid <- unique(c(hosp$ageMid, IFR$ageMid))
    
    ## sample initial values
    logdens <- NA
    while(!is.finite(logdens)) {
        
        eta <- runif(1, 0, 0.01)
        alphaIH <- runif(1, LB, UB)
        alphaID <- runif(1, LB, UB)
        alphaEP <- runif(1, LB, UB)
        alphaHD <- runif(1, LB, UB)
        
        ## set probabilities
        pEP <- exp(alphaEP + eta * ageMid)
        pIH <- exp(alphaIH + eta * ageMid)
        pID <- exp(alphaID + eta * ageMid)
        pHD <- exp(alphaHD + eta * ageMid)
        
        ## validity check
        if(all(pEP <= 1) & all((pIH + pID) <= 1) & all(pHD <= 1)) {
        
            ## log-likelihood estimate using MC sampling
            logdens <- map_dbl(1:nMC, function(j, ageMid, deaths, hosp, pEP, pIH, pID, pHD, nind) {
                
                ## run simulation
                nD <- map(1:length(pEP), function(i, pEP, pIH, pID, pHD, nind) {
                    nP <- rbinom(1, nind, pEP[i])
                    temp <- rmultinom(1, nP, c(pIH[i], pID[i], 1 - pIH[i] - pID[i]))
                    nH <- temp[1]
                    nID <- temp[2]
                    nHD <- rbinom(1, nH, pHD[i])
                    c(nID, nHD, nH)
                }, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
                
                ## extract probabilities
                pD <- map_int(nD, ~sum(.[-3])) / nind
                pH <- map_int(nD, 3) / nind
                
                ## match to data
                hind <- match(hosp$ageMid, ageMid)
                ll <- sum(dbinom(hosp$severe, hosp$all, pH[hind], log = TRUE))
                dind <- match(IFR$ageMid, ageMid)
                ll <- ll + sum(dbinom(IFR$deaths, IFR$all, pD[dind], log = TRUE))
                ll
            }, ageMid = ageMid, deaths = IFR, hosp = hosp, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
            
            ## calculate mean of these replicates on the log-scale
            logdens <- log_sum_exp(logdens)
            
            ## add log-priors
            logdens <- logdens + 
                dunif(alphaIH, LB, UB, log = TRUE) +
                dunif(alphaID, LB, UB, log = TRUE) +
                dunif(alphaEP, LB, UB, log = TRUE) +
                dunif(alphaHD, LB, UB, log = TRUE) +
                dunif(eta, 0, 1, log = TRUE)
            
            ## set up parameters
            pars <- c(alphaEP, alphaID, alphaHD, alphaIH, eta)
        }
    }
    
    ## set up adaptive proposal distribution
    npars <- 5
    propCovIni <- diag(npars) * (0.1^2) / npars
    propCov <- diag(npars) * (0.01^2) / npars
    
    ## set up output
    out <- matrix(NA, N, npars)
    
    ## set up counters
    nacc <- 0
    
    ## run chain
    for(j in 1:N) {
      
        ## propose new values
        u <- runif(1, 0, 1)
        if(u < 0.05) {
            parsProp <- mvrnorm(1, pars, propCovIni)
        } else {
            parsProp <- mvrnorm(1, pars, propCov)
        }
        
        ## transform parameters
        alphaEP <- parsProp[1]
        alphaID <- parsProp[2]
        alphaHD <- parsProp[3]
        alphaIH <- parsProp[4]
        eta <- parsProp[5]
        
        ## set probabilities
        pEP <- exp(alphaEP + eta * ageMid)
        pIH <- exp(alphaIH + eta * ageMid)
        pID <- exp(alphaID + eta * ageMid)
        pHD <- exp(alphaHD + eta * ageMid)
        
        ## add log-priors
        logprop <- dunif(alphaIH, LB, UB, log = TRUE) +
            dunif(alphaID, LB, UB, log = TRUE) +
            dunif(alphaEP, LB, UB, log = TRUE) +
            dunif(alphaHD, LB, UB, log = TRUE) +
            dunif(eta, 0, 1, log = TRUE)
        
        ## validity check
        if(all(pEP <= 1) & all((pIH + pID) <= 1) & all(pHD <= 1) & is.finite(logprop)) {
            
            ## log-likelihood estimate using MC sampling
            temp <- map_dbl(1:nMC, function(j, ageMid, deaths, hosp, pEP, pIH, pID, pHD, nind) {
                
                ## run simulation
                nD <- map(1:length(pEP), function(i, pEP, pIH, pID, pHD, nind) {
                    nP <- rbinom(1, nind, pEP[i])
                    temp <- rmultinom(1, nP, c(pIH[i], pID[i], 1 - pIH[i] - pID[i]))
                    nH <- temp[1]
                    nID <- temp[2]
                    nHD <- rbinom(1, nH, pHD[i])
                    c(nID, nHD, nH)
                }, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
                
                ## extract probabilities
                pD <- map_int(nD, ~sum(.[-3])) / nind
                pH <- map_int(nD, 3) / nind
                
                ## match to data
                hind <- match(hosp$ageMid, ageMid)
                ll <- sum(dbinom(hosp$severe, hosp$all, pH[hind], log = TRUE))
                dind <- match(IFR$ageMid, ageMid)
                ll <- ll + sum(dbinom(IFR$deaths, IFR$all, pD[dind], log = TRUE))
                ll
            }, ageMid = ageMid, deaths = IFR, hosp = hosp, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
            
            ## calculate mean of these replicates on the log-scale
            temp <- log_sum_exp(temp)
            logprop <- logprop + temp
            
            ## accept-reject (symmetric proposal cancels)
            if(log(runif(1, 0, 1)) < (logprop - logdens)) {
                pars <- c(alphaEP, alphaID, alphaHD, alphaIH, eta)
                logdens <- logprop
                nacc <- nacc + 1
            }
        }
        
        ## store samples
        out[j, ] <- pars
        
        ## calculations for adaptive proposal (THIS IS HIGHLY INEFFICIENT AND COULD BE
        ## CALCULATED RECURSIVELY) - uses Sherlock et al. suggestion for pseudo-marginal scaling
        if(j %% nupdate == 0) {
            accrate <- nacc / nupdate
            print(paste0("niter = ", j, ", accrate = ", accrate))
            nacc <- 0
            propCov <- cov(out[min(max(1, j - 25000), N / 4):j, ]) * (2.56^2) / npars
        }
    }
    ## return chain
    colnames(out) <- c("alphaEP", "alphaID", "alphaHD", "alphaIH", "eta")
    out <- as.mcmc(out)
    out
}

## set number of samples and replicates
nsamps <- 100
nreps <- 10

## set seed
set.seed(666)

## set RNG generator to ensure reproducibility
RNGkind("L'Ecuyer-CMRG")

## run MCMCs in parallel
samples <- mclapply(1:(2 * nreps), function(i, hosp, IFR, nsamps) {
  
    ## subsample data
    if(i %% 2 == 0) {
        hosp1 <- filter(hosp, data == "Imperial") %>%
            mutate(all = nsamps) %>%
            mutate(severe = rbinom(rep(1, n()), size = all, prob = propH))
        IFR1 <-  filter(IFR, data == "Imperial") %>%
            mutate(all = nsamps) %>%
            mutate(deaths = rbinom(rep(1, n()), size = all, prob = propD))
    } else {
        hosp1 <- filter(hosp, data == "CDC") %>%
            mutate(all = nsamps) %>%
            mutate(severe = rbinom(rep(1, n()), size = all, prob = propH))
        IFR1 <-  filter(IFR, data == "CDC") %>%
            mutate(all = nsamps) %>%
            mutate(deaths = rbinom(rep(1, n()), size = all, prob = propD))
    }

    ## run MCMC
    samples <- MH(
        N = 1000000, 
        IFR1, 
        hosp1, 
        nind = 100000, 
        nMC = 1, 
        LB = -20, 
        UB = 0
    )
    samples <- window(samples, start = 500000, thin = 10)
    list(samples = samples, hosp1 = hosp1, IFR1 = IFR1)
}, hosp = hosp, IFR = IFR, nsamps = nsamps, mc.cores = 20)

## reset RNG generator
RNGkind("default")

## plot traces (rasterising for space)
map(1:length(samples), function(i, samples) {
    png(paste0("traces", str_pad(i, 2, pad = "0"), "_%d.png"))
    plot(samples[[i]]$samples)
    dev.off()
}, samples = samples)

system("for i in traces*.png; do convert \"$i\" \"$i\".pdf; done")
system("rm traces*.png")
system(paste0("pdftk *.png.pdf output traces.pdf"))
system("rm traces*.png.pdf")

## explore problematic runs
temp <- samples[c(1, 3, 14, 16)]

set.seed(669)

## set RNG generator to ensure reproducibility
RNGkind("L'Ecuyer-CMRG")

## run MCMCs in parallel
temp <- mclapply(1:length(temp), function(i, temp, nsamps) {
    
    ## run MCMC
    samples <- MH(
        N = 1000000, 
        temp[[i]]$IFR1, 
        temp[[i]]$hosp1, 
        nind = 100000, 
        nMC = 1, 
        LB = -20, 
        UB = 0
    )
    samples <- window(samples, start = 500000, thin = 10)
    list(samples = samples, hosp1 = temp[[i]]$hosp1, IFR1 = temp[[i]]$IFR1)
}, temp = temp, nsamps = nsamps, mc.cores = length(temp))

## reset RNG generator
RNGkind("default")

samples[c(1, 3, 14, 16)] <- temp

## plot traces (rasterising for space)
map(1:length(samples), function(i, samples) {
    png(paste0("traces", str_pad(i, 2, pad = "0"), "_%d.png"))
    plot(samples[[i]]$samples)
    dev.off()
}, samples = samples)

system("for i in traces*.png; do convert \"$i\" \"$i\".pdf; done")
system("rm traces*.png")
system(paste0("pdftk *.png.pdf output traces.pdf"))
system("rm traces*.png.pdf")

## union samples
samples <- map(samples, "samples") %>%
    map(as.matrix) %>%
    reduce(rbind)

## ageMids for predictions
ageMid <- sort(unique(c(hosp$ageMid, IFR$ageMid)))

## predictions for deaths
preds <- apply(samples, 1, function(samples, ageMid, nind) {
      
        alphaEP <- samples[1]
        alphaID <- samples[2]
        alphaHD <- samples[3]
        alphaIH <- samples[4]
        eta <- samples[5]
        
        pEP <- exp(alphaEP + eta * ageMid)
        pID <- exp(alphaID + eta * ageMid)
        pHD <- exp(alphaHD + eta * ageMid)
        pIH <- exp(alphaIH + eta * ageMid)
        
        ## check on validity
        if(all(pEP <= 1) & all(pID <= 1) & all(pHD <= 1) & all(pIH <= 1) & all((pID + pIH) <= 1)) {
        
            ## run simulation
            nD <- map(1:length(pEP), function(i, pEP, pIH, pID, pHD, nind) {
                nP <- rbinom(1, nind, pEP[i])
                temp <- rmultinom(1, nP, c(pIH[i], pID[i], 1 - pIH[i] - pID[i]))
                nH <- temp[1]
                nID <- temp[2]
                nHD <- rbinom(1, nH, pHD[i])
                c(nID, nHD, nH)
            }, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
            
            ## extract probabilities
            pD <- map_int(nD, ~sum(.[-3])) / nind
            pH <- map_int(nD, 3) / nind
        } else {
            pD <- rep(NA, length(pEP))
            pH <- rep(NA, length(pEP))
        }
        list(pH, pD)
    }, ageMid = ageMid, nind = 100000)

## extract summaries of predictions
hosp_preds <- map(preds, 1) %>%
    transpose() %>%
    map(unlist) %>%
    map(function(x) {
        x <- x[!is.na(x)]
        tibble(mean = mean(x), LCI = quantile(x, probs = 0.005), UCI = quantile(x, probs = 0.995))
    }) %>%
    bind_rows() %>%
    mutate(ageMid = ageMid) %>%
    mutate(type = "pred") 
death_preds <- map(preds, 2) %>%
    transpose() %>%
    map(unlist) %>%
    map(function(x) {
        x <- x[!is.na(x)]
        tibble(mean = mean(x), LCI = quantile(x, probs = 0.005), UCI = quantile(x, probs = 0.995))
    }) %>%
    bind_rows() %>%
    mutate(ageMid = ageMid) %>%
    mutate(type = "pred") 

## plot predictions
p <- ggplot(death_preds, aes(x = ageMid)) +
    geom_line(aes(y = mean)) +
    geom_point(aes(y = propD), data = IFR) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3)
ggsave("fittedDeaths.pdf", p)
ggsave("fittedDeaths.svg", p)

## plot predictions
p <- ggplot(hosp_preds, aes(x = ageMid)) +
    geom_line(aes(y = mean)) +
    geom_point(aes(y = propH), data = hosp) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3)
ggsave("fittedHosp.pdf", p)
ggsave("fittedHosp.svg", p)

## fit range of finite mixture models
samples1 <- samples[sample.int(nrow(samples), 50000), ]
# samples1 <- maximin(samples, 10)
mod <- densityMclust(samples1)

## summary of finite mixture models
summary(mod)
plot(mod, what = "BIC")

## take random samples from mixture
nimp <- 50000
props <- sim(mod$modelName, mod$parameters, nimp)[, -1] %>%
    as_tibble() %>%
    set_names(c("alphaEP", "alphaID", "alphaHD", "alphaIH", "eta")) %>%
    mutate(Estimate = "Mixture") %>%
    filter(eta > 0)

## plot mixture posterior against initial posterior
p <- as_tibble(samples1) %>%
    mutate(Estimate = "MCMC") %>%
    rbind(props) %>%
    ggpairs(mapping = aes(colour = Estimate, alpha = 0.5), upper = list(continuous = "density"), columns = 1:5, legend = 1) + 
    theme(legend.position = "bottom")
leg <- getPlot(p, 1, 1) +
    guides(alpha = FALSE)
leg <- grab_legend(leg)
p <- as_tibble(samples1) %>%
    mutate(Estimate = "MCMC") %>%
    rbind(props) %>%
    ggpairs(mapping = aes(colour = Estimate, alpha = 0.5), upper = list(continuous = "density"), columns = 1:5, legend = leg) + 
    theme(legend.position = "bottom")
ggsave("simsMixturePosterior.png", p)

## save mclust object
saveRDS(mod, "fullmod.rds")
