## load libraries
library(mclust)
library(MASS)
library(tidyverse)
library(coda)
library(GGally)
library(svglite)

## read in reference estimates
IFR <- readRDS("IFR.rds")

## transform according to bounds
IFR <- mutate(IFR, all = round(IFR * (1 - IFR) * (1.96 / ((UCI - LCI) / 2))^2)) %>%
  mutate(severe = round(IFR * all))

## read in FMM
model <- readRDS("hospitalisation.rds")

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
MH <- function(N, model, ageMid, severe, all, nind, nMC, LB = -20) {
    
    ## sample initial values
    logdens <- NA
    while(!is.finite(logdens)) {
        eta <- sim(model$modelName, model$parameters, 1)[, -1]
        alphaIH <- eta[1]
        eta <- eta[2]
        pIH <- exp(alphaIH + eta * ageMid)
        stopifnot(all(pIH <= 1))
        
        alphaID <- runif(1, LB, log(1 - exp(alphaIH + eta * max(ageMid))) - eta * max(ageMid))
        pID <- exp(alphaID + eta * ageMid)
        stopifnot(all(pID <= 1))
        
        alphaEP <- runif(1, LB, -eta * max(ageMid))
        pEP <- exp(alphaEP + eta * ageMid)
        stopifnot(all(pEP <= 1))
        
        alphaHD <- runif(1, LB, -eta * max(ageMid))
        pHD <- exp(alphaHD + eta * ageMid)
        stopifnot(all(pHD <= 1))
        
        ## log-likelihood estimate using MC sampling
        logdens <- map_dbl(1:nMC, function(j, ageMid, pEP, pIH, pID, pHD, nind, severe, all) {
            nD <- map_int(1:length(ageMid), function(i, pEP, pIH, pID, pHD, nind) {
              nP <- rbinom(1, nind, pEP[i])
              temp <- rmultinom(1, nP, c(pIH[i], pID[i], 1 - pIH[i] - pID[i]))
              nH <- temp[1]
              nID <- temp[2]
              nHD <- rbinom(1, nH, pHD[i])
              nID + nHD
          }, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
          pD <- nD / nind
          sum(dbinom(severe, size = all, prob = pD, log = TRUE))
        }, ageMid = ageMid, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind, severe = severe, all = all)
        
        ## calculate mean of these replicates on the log-scale
        logdens <- log_sum_exp(logdens)
        
        ## add log-priors
        logdens <- logdens + 
          dens(model$modelName, matrix(c(alphaIH, eta), nrow = 1), model$parameters, logarithm = TRUE) +
          dunif(alphaID, LB, log(1 - exp(alphaIH + eta * max(ageMid))) - eta * max(ageMid), log = TRUE) +
          dunif(alphaEP, LB, -eta * max(ageMid), log = TRUE) +
          dunif(alphaHD, LB, -eta * max(ageMid), log = TRUE)
        
        ## check eta (don't need explicit truncation normalisation because it cancels)
        if(eta < 0) logdens <- NA
        
        ## set up parameters
        pars <- c(alphaEP, alphaID, alphaHD, alphaIH, eta)
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
        pEP <- exp(alphaEP + eta * ageMid)
        pID <- exp(alphaID + eta * ageMid)
        pHD <- exp(alphaHD + eta * ageMid)
        pIH <- exp(alphaIH + eta * ageMid)
        
        ## check log-priors
        logprop <- dens(model$modelName, matrix(c(alphaIH, eta), nrow = 1), model$parameters, logarithm = TRUE) +
          dunif(alphaID, LB, log(1 - exp(alphaIH + eta * max(ageMid))) - eta * max(ageMid), log = TRUE) +
          dunif(alphaEP, LB, -eta * max(ageMid), log = TRUE) +
          dunif(alphaHD, LB, -eta * max(ageMid), log = TRUE)
        
        ## check eta (don't need explicit truncation normalisation because it cancels)
        if(eta < 0) logprop <- NA
    
        ## log-likelihood estimate
        if(is.finite(logprop)) {
            ## quick additional check
            stopifnot(all(pEP <= 1) & all(pID <= 1) & all(pHD <= 1) & all(pIH <= 1) & all((pID + pIH) <= 1))
          
            ## log-likelihood estimate using MC sampling
            temp <- map_dbl(1:nMC, function(j, ageMid, pEP, pIH, pID, pHD, nind, severe, all) {
              nD <- map_int(1:length(ageMid), function(i, pEP, pIH, pID, pHD, nind) {
                nP <- rbinom(1, nind, pEP[i])
                temp <- rmultinom(1, nP, c(pIH[i], pID[i], 1 - pIH[i] - pID[i]))
                nH <- temp[1]
                nID <- temp[2]
                nHD <- rbinom(1, nH, pHD[i])
                nID + nHD
              }, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
              pD <- nD / nind
              sum(dbinom(severe, size = all, prob = pD, log = TRUE))
            }, ageMid = ageMid, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind, severe = severe, all = all)
            
            ## calculate mean of these replicates
            temp <- log_sum_exp(temp)
            logprop <- logprop + temp
            
            ## accept-reject (symmetric proposal cancels)
            if(log(runif(1, 0, 1)) < (logprop - logdens)) {
                pars <- c(alphaEP, alphaID, alphaHD, alphaIH, eta)
                logdens <- logprop
                nacc <- nacc + 1
            }
        }
        
        # store samples
        out[j, ] <- pars
        
        ## calculations for adaptive proposal (THIS IS HIGHLY INEFFICIENT AND COULD BE
        ## CALCULATED RECURSIVELY) - uses Sherlock et al. suggestion for pseudo-marginal
        if(j %% 100 == 0) {
            accrate <- nacc / 100
            print(paste0("j = ", j, ", accrate = ", accrate))
            nacc <- 0
            propCov <- cov(out[1:j, ]) * (2.56^2) / npars
        }
    }
    ## return chain
    colnames(out) <- c("alphaEP", "alphaID", "alphaHD", "alphaIH", "eta")
    out <- as.mcmc(out)
    out
}

## run MCMC
samples <- MH(50000, model, IFR$ageMid, IFR$severe, IFR$all, 10000, 10)

## plot traces
pdf("simtraces.pdf")
plot(samples)
dev.off()

## plot posterior
p <- as.matrix(window(samples, start = 25000)) %>%
  as_tibble() %>%
  ggpairs(upper = list(continuous = "density"))
ggsave("simsPosterior.png", p)

## predictions
preds <- apply(as.matrix(window(samples, start = 25000)), 1, function(samples, ageMid, nind) {
  
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
    stopifnot(all(pEP <= 1) & all(pID <= 1) & all(pHD <= 1) & all(pIH <= 1) & all((pID + pIH) <= 1))
    
    ## run simulation
    nD <- map_int(1:length(ageMid), function(i, pEP, pIH, pID, pHD, nind) {
        nP <- rbinom(1, nind, pEP[i])
        temp <- rmultinom(1, nP, c(pIH[i], pID[i], 1 - pIH[i] - pID[i]))
        nH <- temp[1]
        nID <- temp[2]
        nHD <- rbinom(1, nH, pHD[i])
        nID + nHD
      }, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
    nD / nind
}, ageMid = IFR$ageMid, nind = 100000)
preds <- apply(preds, 1, function(x) {
        tibble(mean = mean(x), LCI = quantile(x, probs = 0.005), UCI = quantile(x, probs = 0.995))
    }) %>%
  bind_rows() %>%
  mutate(ageMid = IFR$ageMid)%>%
  mutate(type = "pred") 

## plot predictions
p <- IFR %>%
    select(mean = IFR, LCI, UCI, ageMid) %>%
    mutate(type = "IFR") %>%
    rbind(preds) %>%
    ggplot(aes(x = ageMid)) +
      geom_point(aes(y = mean, colour = type)) +
      geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = type), alpha = 0.3)
ggsave("simsPreds.pdf", p)
ggsave("simsPreds.svg", p)
ggsave("simsPreds.png", p)
