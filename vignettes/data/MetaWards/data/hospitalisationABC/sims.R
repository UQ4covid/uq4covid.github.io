## load libraries
library(mclust)
library(MASS)
library(tidyverse)
library(svglite)
library(SimBIID)

## read in data
IFR <- readRDS("IFR.rds")
hosp <- readRDS("hosp.rds")

## generate function to run simulator
## and return summary statistics
simABC <- function(pars, data, tols, u, u1) {
  
    u <- u1
    
    ## extract parameters
    alphaEP <- pars[1]
    alphaID <- pars[2]
    alphaHD <- pars[3]
    alphaIH <- pars[4]
    eta <- pars[5]
    
    ## rudimentary check of inputs
    stopifnot(names(u) %in% c("nind", "ageMidH", "ageMidD", "ageMax"))
    
    ## set max probabilities
    pEP <- exp(alphaEP + eta * u$ageMax)
    pIH <- exp(alphaIH + eta * u$ageMax)
    pID <- exp(alphaID + eta * u$ageMax)
    pHD <- exp(alphaHD + eta * u$ageMax)
    
    ## reject if invalid (since simulation will reject anyway)
    if(pEP > 1 | pIH > 1 | pHD > 1 | (pID + pIH) > 1) {
        return(NA)
    }
    
    ## set probabilities for age classes
    pIH <- exp(alphaIH + eta * u$ageMidH)
    
    ## simulate hospitalisation data
    pH <- rbinom(rep(1, length(pIH)), size = rep(u$nind, length(pIH)), pIH) / u$nind
    
    ## set probabilities for age classes
    pEP <- exp(alphaEP + eta * u$ageMidD)
    pIH <- exp(alphaIH + eta * u$ageMidD)
    pID <- exp(alphaID + eta * u$ageMidD)
    pHD <- exp(alphaHD + eta * u$ageMidD)
  
    ## simulate from latent variable model for each age-class
    nD <- map_int(1:length(pEP), function(i, pEP, pIH, pID, pHD, nind) {
      nP <- rbinom(1, nind, pEP[i])
      temp <- rmultinom(1, nP, c(pIH[i], pID[i], 1 - pIH[i] - pID[i]))
      nH <- temp[1]
      nID <- temp[2]
      nHD <- rbinom(1, nH, pHD[i])
      nID + nHD
    }, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = u$nind)
    pD <- nD / u$nind
    
    ## return vector if match, else return NA
    if(all(abs(c(pH, pD) - data) <= tols)) {
        return(c(pH, pD))
    } else {
        return(NA)
    }
}

## set priors
priors <- data.frame(
  parnames = c("alphaEP", "alphaID", "alphaIH", "alphaHD", "eta"), 
  dist = rep("unif", 5), 
  stringsAsFactors = FALSE
)
priors$p1 <- c(-20, -20, -20, -20, 0)
priors$p2 <- c(0, 0, 0, 0, 0.1)

## define the targeted summary statistics
data <- c(hosp$propH, IFR$propD)
names(data) <- c(paste0("H", 1:length(hosp$ageMid)), paste0("D", 1:length(IFR$ageMid)))

## set other inputs
u1 <- list(nind = 100000, ageMidH = hosp$ageMid, ageMidD = IFR $ageMid, ageMax = 100)

## set initial tolerances
tols <- rep(0.2, length(data))
colnames(tols) <- names(data)

## run generations of ABC-SMC
post <- ABCSMC(
  x = data, 
  priors = priors, 
  func = simABC, 
  u = 1:3, 
  tols = tols[1, ],
  ptols = 0.5,
  ngen = 3,
  npart = 5000,
  u1 = u1,
  parallel = TRUE
)

plot(post)
plot(post, "output")

## run for some more generations of ABC
post <- ABCSMC(
  x = post,
  ptols = 0.5,
  ngen = 5,
  parallel = TRUE
)
plot(post, "output")

## run for some more generations of ABC
post <- ABCSMC(
  x = post,
  ptols = 0.2,
  ngen = 5,
  parallel = TRUE
)
plot(post, "output", gen = 13)




## M-H algorithm
MH <- function(N, IFR, hosp, nind, nMC, LB = -20, maxAge = 100, nupdate = 100) {
  
    ## split data
    hospCDC <- filter(hosp, data == "CDC") %>%
        select(hosp, all, ageMid)
    hospImp <- filter(hosp, data == "Imperial") %>%
        select(hosp, all, ageMid)
    
    ## sample initial values
    logdens <- NA
    while(!is.finite(logdens)) {
        eta <- runif(1, 0, 0.01)
        alphaIH <- runif(1, LB, -eta * maxAge)
        alphaID <- runif(1, LB, log(1 - exp(alphaIH + eta * maxAge)) - eta * maxAge)
        alphaEP <- runif(1, LB, -eta * maxAge)
        alphaHD <- runif(1, LB, -eta * maxAge)
        
        ## log-likelihood for hospitalisation data x 2
        logdens <- sum(dbinom(hospImp$hosp, hospImp$all, exp(alphaIH + eta * hospImp$ageMid), log = TRUE))
        logdens <- logdens + sum(dbinom(hospCDC$hosp, hospCDC$all, exp(alphaIH + eta * hospCDC$ageMid), log = TRUE))
        
        ## set probabilities
        pEP <- exp(alphaEP + eta * IFR$ageMid)
        pIH <- exp(alphaIH + eta * IFR$ageMid)
        pID <- exp(alphaID + eta * IFR$ageMid)
        pHD <- exp(alphaHD + eta * IFR$ageMid)
        
        ## log-likelihood estimate using MC sampling
        temp <- map_dbl(1:nMC, function(j, ageMid, deaths, all, pEP, pIH, pID, pHD, nind) {
            nD <- map_int(1:length(pEP), function(i, pEP, pIH, pID, pHD, nind) {
              nP <- rbinom(1, nind, pEP[i])
              temp <- rmultinom(1, nP, c(pIH[i], pID[i], 1 - pIH[i] - pID[i]))
              nH <- temp[1]
              nID <- temp[2]
              nHD <- rbinom(1, nH, pHD[i])
              nID + nHD
          }, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
          pD <- nD / nind
          sum(dbinom(deaths, all, pD, log = TRUE))
        }, ageMid = IFR$ageMid, deaths = IFR$deaths, all = IFR$all, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
        
        ## calculate mean of these replicates on the log-scale
        temp <- log_sum_exp(temp)
        
        ## add log-priors
        logdens <- logdens + temp + 
          dunif(alphaIH, LB, -eta * maxAge, log = TRUE) +
          dunif(alphaID, LB, log(1 - exp(alphaIH + eta * maxAge)) - eta * maxAge, log = TRUE) +
          dunif(alphaEP, LB, -eta * maxAge, log = TRUE) +
          dunif(alphaHD, LB, -eta * maxAge, log = TRUE) +
          dunif(eta, 0, 1, log = TRUE)
        
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
        
        ## check valid according to priors
        logprop <- dunif(alphaIH, LB, -eta * maxAge, log = TRUE) +
            dunif(alphaID, LB, log(1 - exp(alphaIH + eta * maxAge)) - eta * maxAge, log = TRUE) +
            dunif(alphaEP, LB, -eta * maxAge, log = TRUE) +
            dunif(alphaHD, LB, -eta * maxAge, log = TRUE) +
            dunif(eta, 0, 1, log = TRUE)
        
        ## log-likelihood estimate
        if(is.finite(logprop)) {
          
            ## log-likelihood for hospitalisation data x 2
            logprop <- logprop + sum(dbinom(hospImp$hosp, hospImp$all, exp(alphaIH + eta * hospImp$ageMid), log = TRUE))
            logprop <- logprop + sum(dbinom(hospCDC$hosp, hospCDC$all, exp(alphaIH + eta * hospCDC$ageMid), log = TRUE))
            
            ## set probabilities
            pEP <- exp(alphaEP + eta * IFR$ageMid)
            pIH <- exp(alphaIH + eta * IFR$ageMid)
            pID <- exp(alphaID + eta * IFR$ageMid)
            pHD <- exp(alphaHD + eta * IFR$ageMid)
            
            ## log-likelihood estimate using MC sampling
            temp <- map_dbl(1:nMC, function(j, ageMid, deaths, all, pEP, pIH, pID, pHD, nind) {
              nD <- map_int(1:length(pEP), function(i, pEP, pIH, pID, pHD, nind) {
                nP <- rbinom(1, nind, pEP[i])
                temp <- rmultinom(1, nP, c(pIH[i], pID[i], 1 - pIH[i] - pID[i]))
                nH <- temp[1]
                nID <- temp[2]
                nHD <- rbinom(1, nH, pHD[i])
                nID + nHD
              }, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
              pD <- nD / nind
              sum(dbinom(deaths, all, pD, log = TRUE))
            }, ageMid = IFR$ageMid, deaths = IFR$deaths, all = IFR$all, pEP = pEP, pIH = pIH, pID = pID, pHD = pHD, nind = nind)
            
            ## calculate mean of these replicates on the log-scale
            temp <- log_sum_exp(temp)
            
            ## add log-priors
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
            propCov <- cov(out[min(max(1, j - 10000), N / 2):j, ]) * (2.56^2) / npars
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
samples <- list()
pdf("traces.pdf")
for(i in 1:nreps) {
  
    ## subsample data
    hosp1 <- mutate(hosp, all = nsamps) %>%
       mutate(hosp = rbinom(rep(1, n()), size = all, prob = propH))
    IFR1 <- mutate(IFR, all = nsamps) %>%
        mutate(deaths = rbinom(rep(1, n()), size = all, prob = propD))

    ## run MCMC
    temp <- MH(
        N = 200000, 
        IFR1, 
        hosp1, 
        nind = 100000, 
        nMC = 1, 
        LB = -20, 
        maxAge = 100,
        nupdate = 1000
    )
    temp <- window(temp, start = 100000, thin = 4)
    plot(temp)
    samples[[i]] <- as.matrix(temp)
}
dev.off()

## union samples
samples <- reduce(samples, rbind)

## hospitalisation predictions
hosp_preds <- matrix(1, nrow = 1, ncol = nrow(hosp)) %>%
    rbind(hosp$ageMid)
hosp_preds <- exp(samples[, match(c("alphaIH", "eta"), colnames(samples))] %*% hosp_preds)
hosp_preds1 <- rbind(rep(nsamps, ncol(hosp_preds)), hosp_preds) %>%
    apply(2, function(x) {
      rbinom(rep(1, length(x) - 1), x[1], x[-1]) / x[1]
    })

## summarise predictions
hosp_preds_sum <- hosp_preds %>%
  apply(2, function(x) {
    tibble(mean = mean(x), LCI = quantile(x, probs = 0.005), UCI = quantile(x, probs = 0.995))
  }) %>%
  bind_rows() %>%
  mutate(ageMid = hosp$ageMid) %>%
  arrange(ageMid)
hosp_preds1_sum <- hosp_preds1 %>%
  apply(2, function(x) {
    tibble(mean = mean(x), LCI = quantile(x, probs = 0.005), UCI = quantile(x, probs = 0.995))
  }) %>%
  bind_rows() %>%
  mutate(ageMid = hosp$ageMid) %>%
  arrange(ageMid)

## fitted plot
p <- ggplot(hosp, aes(x = ageMid)) +
  geom_point(aes(y = propH)) +
  geom_line(aes(y = mean), data = hosp_preds_sum) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), data = hosp_preds_sum, alpha = 0.5) +
  geom_ribbon(aes(ymin = LCI, ymax = UCI), data = hosp_preds1_sum, alpha = 0.5) +
  xlab("Age") + ylab("Probability of hospitalisation")
ggsave("fittedhosp.pdf", p)
ggsave("fittedhosp.svg", p)

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
p <- ggplot(preds, aes(x = ageMid)) +
      geom_line(aes(y = mean)) +
      geom_point(aes(y = propD), data = IFR) +
      geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.3)
ggsave("fittedDeaths.pdf", p)
ggsave("fittedDeaths.svg", p)
ggsave("fittedDeaths.png", p)

## fit range of finite mixture models
samples <- samples[sample.int(nrow(samples), 30000), ]
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

## save mclust object
saveRDS(mod, "fullmod.rds")
