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
IFR <- mutate(IFR, diffL = IFR - LCI) %>%
    mutate(diffU = UCI - IFR) %>%
    mutate(diff = ifelse(diffL > diffU, diffL, diffU)) %>%
    mutate(diff = (diff / 1.96)^2) %>%
    mutate(all = round(IFR * (1 - IFR) / diff)) %>%
    mutate(severe = round(IFR * all)) %>%
    select(-starts_with("diff"))

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
    propCov <- diag(npars) * (0.1^2) / npars
    
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
        
        ## store samples
        out[j, ] <- pars
        
        ## calculations for adaptive proposal (THIS IS HIGHLY INEFFICIENT AND COULD BE
        ## CALCULATED RECURSIVELY) - uses Sherlock et al. suggestion for pseudo-marginal scaling
        if(j %% 100 == 0) {
            accrate <- nacc / 100
            print(paste0("niter = ", j, ", accrate = ", accrate))
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
samples <- MH(100000, model, IFR$ageMid, IFR$severe, IFR$all, 10000, 10)
samples <- window(samples, start = 10000, thin = 3)

## plot traces
pdf("simtraces.pdf")
plot(samples)
dev.off()

## convert to matrix
samples <- as.matrix(samples)

## predictions
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
    }, ageMid = IFR$ageMid, nind = 100000) %>%
    apply(1, function(x) {
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

## fit range of finite mixture models
mod <- densityMclust(samples)

## summary of finite mixture models
summary(mod)
plot(mod, what = "BIC")

## take random samples from mixture
nimp <- 30000
props <- sim(mod$modelName, mod$parameters, nimp)[, -1] %>%
    as_tibble() %>%
    set_names(c("alphaEP", "alphaID", "alphaHD", "alphaIH", "eta")) %>%
    mutate(Estimate = "Mixture") %>%
    filter(eta > 0)

## plot mixture posterior against initial posterior
p <- as_tibble(samples) %>%
    mutate(Estimate = "MCMC") %>%
    rbind(props) %>%
    ggpairs(mapping = aes(colour = Estimate, alpha = 0.5), upper = list(continuous = "density"), columns = 1:5, legend = 1) + 
    theme(legend.position = "bottom")
leg <- getPlot(p, 1, 1) +
    guides(alpha = FALSE)
leg <- grab_legend(leg)
p <- as_tibble(samples) %>%
    mutate(Estimate = "MCMC") %>%
    rbind(props) %>%
    ggpairs(mapping = aes(colour = Estimate, alpha = 0.5), upper = list(continuous = "density"), columns = 1:5, legend = leg) + 
    theme(legend.position = "bottom")
ggsave("simsMixturePosterior.png", p)
