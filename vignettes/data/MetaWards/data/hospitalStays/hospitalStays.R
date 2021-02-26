## load libraries
library(mclust)
library(tidyverse)
library(fitdistrplus)
library(nimble)
library(coda)
library(GGally)
library(svglite)

## load data
hosp <- read_csv("admissionToOutcomeLineList.csv", col_names = TRUE)

## summarise data
mutate_if(hosp, is.character, as.factor) %>%
    summary()

## match data roughly to age-classes used in the model
hosp <- mutate(hosp, ageCat = ifelse(ageCat == "70-79", "70+", ageCat)) %>%
    mutate(ageCat = ifelse(ageCat == "80-89", "70+", ageCat)) %>%
    mutate(ageCat = ifelse(ageCat == "90+", "70+", ageCat)) %>%
    arrange(ageCat)

## fit some simple models to each age-class and ignore outcome
hosp <- hosp %>%
    filter(finaloutcome == "Death" | finaloutcome == "Discharged") %>%
    dplyr::select(time, ageCat) %>%
    group_by(ageCat) %>%
    nest() %>%
    mutate(exp = map(data, ~fitdist(.$time, "exp"))) %>%
    mutate(gamma = map(data, ~fitdist(.$time, "gamma"))) %>%
    mutate(expAIC = map_dbl(exp, "aic")) %>%
    mutate(gammaAIC = map_dbl(gamma, "aic")) %>%
    mutate(model = ifelse(expAIC < gammaAIC, exp, gamma)) %>%
    mutate(coefs = map(model, "estimate")) %>%
    mutate(preds = map(coefs, ~{
        t <- 0:100
        if(length(.) == 1) {
            out <- dexp(t, rate = .)
        } else {
            out <- dgamma(t, shape = .[1], rate = .[2])
        }
        tibble(time = t, out = out)
    }))

## produce predictions using exponential model
hosp <- hosp %>%
    mutate(coefs = map_dbl(exp, "estimate")) %>%
    mutate(preds = map(coefs, ~{
        t <- 0:100
        out <- dexp(t, rate = .)
        tibble(time = t, out = out)
    }))

## fitted plots
p <- dplyr::select(hosp, ageCat, data) %>%
    unnest(cols = data) %>%
    ggplot(aes(x = time)) +
    geom_density() +
    facet_wrap(~ageCat, scales = "free_y") +
    geom_line(aes(y = out), data = dplyr::select(hosp, ageCat, preds) %>%
                  unnest(cols = preds), col = "red") +
    xlab("Time") + ylab("Density")
ggsave("expFits.pdf", p)
ggsave("expFits.svg", p)

## extract mid-points of age-classes
hosp <- hosp %>%
    mutate(ageCat = ifelse(ageCat == "70+", "70-80", ageCat)) %>%
    separate(ageCat, c("LB", "UB"), sep = "-", remove = FALSE) %>%
    mutate_at(vars(c("UB", "LB")), as.numeric) %>%
    mutate(ageMid = (LB + UB) / 2) %>%
    dplyr::select(-LB, -UB) %>%
    mutate(TH = 1 /coefs)

## fit Bayesian hierarchical model to account for sampling variability
code <- nimbleCode({
    
    ## fit exponential survival model to individual data
    ## stratified by age
    for (i in 1:nInd) {
        t[i] ~ dexp(lambda[ageGrp[i]])
    }
    
    ## fit log-linear model to 1 / rates
    for(k in 1:nAge) {
        lambda[k] <- 1.0 / TH[k]
        mu[k] <- alpha + eta * ageMid[k]
        log(TH[k]) ~ dnorm(mu[k], sigma[k])
        log(sigma[k]) ~ dnorm(muS, sigmaS)
    }
    
    ## priors
    alpha ~ dnorm(0, 5)
    eta ~ T(dnorm(0, 1), 0, )
    muS ~ dnorm(0, 5)
    sigmaS ~ dexp(1)
})

## set up other components of model
hospBayes <- hosp %>%
    dplyr::select(ageCat, data, ageMid) %>%
    unnest(data) %>%
    mutate(ageGrp = factor(ageCat)) %>%
    mutate(ageGrp = as.numeric(ageGrp))
consts <- list(
    nInd = nrow(hospBayes), 
    nAge = length(unique(hospBayes$ageGrp)),
    ageGrp = hospBayes$ageGrp,
    ageMid = sort(unique(hospBayes$ageMid)),
    maxAge = 100)
data <- list(t = hospBayes$time)

## sample initial values
initFn <- function(ageMid) {
    alpha <- rnorm(1, 0, 1)
    eta <- abs(rnorm(1, 0, 1))
    sigma <- rexp(length(ageMid), 1)
    sigmaS <- rexp(1, 1)
    muS <- rnorm(1, 0, 1)
    log_TH <- rnorm(alpha + eta * ageMid, sigma)
    TH <- exp(log_TH)
    lambda <- 1 / TH
    inits <- list(
        alpha = alpha,
        eta = eta,
        sigma = sigma,
        lambda = lambda,
        sigmaS = sigmaS,
        muS = muS,
        TH = TH,
        log_TH = log_TH
    )
    inits
}

## define the model, data, inits and constants
model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(consts$ageMid))

## compile the model
cmodel <- compileNimble(model)

## set monitors
config <- configureMCMC(cmodel, monitors = c("alpha", "eta", "lambda"), thin = 1)

## check monitors and samplers
config$printMonitors()
config$printSamplers()

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run <- runMCMC(cbuilt, 
    niter = 20000, 
    nburnin = 5000, 
    nchains = 2, 
    progressBar = TRUE, 
    summary = TRUE, 
    samplesAsCodaMCMC = TRUE, 
    thin = 1))

## plot traces
pdf("traces.pdf")
plot(run$samples)
dev.off()
samples <- as.matrix(run$samples)

## predictions
hosp_preds <- matrix(1, nrow = 1, ncol = 100) %>%
    rbind(seq(min(hosp$ageMid), max(hosp$ageMid), length.out = 100))
hosp_preds <- samples[, 1:2] %*% hosp_preds
hosp_preds <- apply(hosp_preds, 2, function(x) {
        tibble(mean = mean(x), LCI = quantile(x, probs = 0.005), UCI = quantile(x, probs = 0.995))
    }) %>%
    bind_rows() %>%
    mutate(ageMid = seq(min(hosp$ageMid), max(hosp$ageMid), length.out = 100))

## fitted plot
p <- ggplot(hosp, aes(x = ageMid)) +
    geom_point(aes(y = log(TH))) +
    geom_line(aes(y = mean), data = hosp_preds) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), data = hosp_preds, alpha = 0.5) +
    xlab("Age") + ylab("log(mean hospital stay time)")# +
    # geom_abline(intercept = 1.36234700, slope =  0.01560438, colour = "red")
ggsave("fittedLnBayes.pdf", p)
ggsave("fittedLnBayes.svg", p)

## posterior correlation
cor(samples[, 1:2])

## fit range of finite mixture models
mod <- densityMclust(samples[, 1:2])

## summary of finite mixture models
summary(mod)
plot(mod, what = "BIC")

## take random samples from mixture
nimp <- 30000
props <- sim(mod$modelName, mod$parameters, nimp)[, -1] %>%
    as_tibble() %>%
    set_names(c("alpha", "eta")) %>%
    mutate(Estimate = "Mixture") %>%
    filter(eta > 0)

p <- as_tibble(samples[, 1:2]) %>%
    mutate(Estimate = "MCMC") %>%
    rbind(props) %>%
    ggpairs(mapping = aes(colour = Estimate, alpha = 0.5), upper = list(continuous = "density"), columns = 1:2, legend = 1) + 
    theme(legend.position = "bottom")
leg <- getPlot(p, 1, 1) +
    guides(alpha = FALSE)
leg <- grab_legend(leg)
p <- as_tibble(samples[, 1:2]) %>%
    mutate(Estimate = "MCMC") %>%
    rbind(props) %>%
    ggpairs(mapping = aes(colour = Estimate, alpha = 0.5), upper = list(continuous = "density"), columns = 1:2, legend = leg) + 
    theme(legend.position = "bottom")
ggsave("mixturePosterior.pdf", p)
ggsave("mixturePosterior.svg", p)
ggsave("mixturePosterior.png", p)

## save mclust object
saveRDS(mod, "hospStays.rds")
