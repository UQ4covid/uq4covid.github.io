## load libraries
library(mclust)
library(tidyverse)
library(nimble)
library(coda)
library(GGally)
library(svglite)
library(readxl)

## load data
hosp <- readRDS("hosp.rds")

## fit Bayesian model
code <- nimbleCode({
    
    ## fit log-linear model to 1 / rates
    for(k in 1:nAge) {
        severe[k] ~ dbin(pS[k], all[k])
        log(pS[k]) <- alpha + eta * ageMid[k]
    }
    
    ## priors
    alpha ~ T(dnorm(0, sd = 5), , -eta * 100)
    eta ~ T(dnorm(0, sd = 1), 0, )
})

## sample initial values
initFn <- function(ageMid) {
    pS <- 2
    while(any(pS >= 1 | pS <= 0)) {
        alpha <- rnorm(1, 0, 1)
        eta <- abs(rnorm(1, 0, 1))
        pS <- exp(alpha + eta * ageMid)
    }
    inits <- list(
        alpha = alpha,
        eta = eta,
        pS = pS
    )
    inits
}

nrep <- 10
nsamps <- 100
j <- 1
samples <- list()
pdf("traces.pdf")
for(i in 1:nrep) {
    
    ## subsample CDC
    hosp1 <- filter(hosp, data == "CDC") %>%
        mutate(all = nsamps) %>%
        mutate(severe = rbinom(rep(1, n()), size = all, prob = prop))

    ## set up other components of model
    consts <- list(
        nAge = nrow(hosp1),
        ageMid = hosp1$ageMid,
        all = hosp1$all)
    data <- list(severe = hosp1$severe)
    
    ## define the model, data, inits and constants
    model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(100))
    
    ## compile the model
    cmodel <- compileNimble(model)
    
    ## set monitors
    config <- configureMCMC(cmodel, monitors = c("alpha", "eta"), thin = 1)
    
    ## check monitors and samplers
    config$printMonitors()
    config$printSamplers()
    
    ## build the model
    built <- buildMCMC(config)
    cbuilt <- compileNimble(built)
    
    ## run the model
    system.time(run <- runMCMC(cbuilt, 
        niter = 155000, 
        nburnin = 5000, 
        nchains = 2, 
        progressBar = TRUE, 
        summary = TRUE, 
        samplesAsCodaMCMC = TRUE, 
        thin = 10))

    ## plot traces
    plot(run$samples)
    
    ## extract samples
    samples[[j]] <- as.matrix(run$samples)
    j <- j + 1
    
    ## subsample CDC and remove final age class
    hosp1 <- filter(hosp, data == "CDC") %>%
        filter(ageMid != max(ageMid)) %>%
        mutate(all = nsamps) %>%
        mutate(severe = rbinom(rep(1, n()), size = all, prob = prop))
    
    ## set up other components of model
    consts <- list(
        nAge = nrow(hosp1),
        ageMid = hosp1$ageMid,
        all = hosp1$all)
    data <- list(severe = hosp1$severe)
    
    ## define the model, data, inits and constants
    model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(consts$ageMid))
    
    ## compile the model
    cmodel <- compileNimble(model)
    
    ## set monitors
    config <- configureMCMC(cmodel, monitors = c("alpha", "eta"), thin = 1)
    
    ## check monitors and samplers
    config$printMonitors()
    config$printSamplers()
    
    ## build the model
    built <- buildMCMC(config)
    cbuilt <- compileNimble(built)
    
    ## run the model
    system.time(run <- runMCMC(cbuilt, 
        niter = 155000, 
        nburnin = 5000, 
        nchains = 2, 
        progressBar = TRUE, 
        summary = TRUE, 
        samplesAsCodaMCMC = TRUE, 
        thin = 10))
    
    ## plot traces
    plot(run$samples)
    
    ## extract samples
    samples[[j]] <- as.matrix(run$samples)
    j <- j + 1
    
    ## subsample Imperial
    hosp1 <- filter(hosp, data == "Imperial") %>%
        mutate(all = nsamps) %>%
        mutate(severe = rbinom(rep(1, n()), size = all, prob = prop))
    
    ## set up other components of model
    consts <- list(
        nAge = nrow(hosp1),
        ageMid = hosp1$ageMid,
        all = hosp1$all)
    data <- list(severe = hosp1$severe)
    
    ## define the model, data, inits and constants
    model <- nimbleModel(code = code, constants = consts, data = data, inits = initFn(consts$ageMid))
    
    ## compile the model
    cmodel <- compileNimble(model)
    
    ## set monitors
    config <- configureMCMC(cmodel, monitors = c("alpha", "eta"), thin = 1)
    
    ## check monitors and samplers
    config$printMonitors()
    config$printSamplers()
    
    ## build the model
    built <- buildMCMC(config)
    cbuilt <- compileNimble(built)
    
    ## run the model
    system.time(run <- runMCMC(cbuilt, 
        niter = 155000, 
        nburnin = 5000, 
        nchains = 2, 
        progressBar = TRUE, 
        summary = TRUE, 
        samplesAsCodaMCMC = TRUE, 
        thin = 10))
    
    ## plot traces
    plot(run$samples)
    
    ## extract samples
    samples[[j]] <- as.matrix(run$samples)
    j <- j + 1
}
dev.off()

## union samples
samples <- reduce(samples, rbind)

## order hosp
hosp1 <- arrange(hosp, ageMid) %>%
    distinct(ageMid, .keep_all = TRUE)

## predictions
hosp_preds <- matrix(1, nrow = 1, ncol = nrow(hosp1)) %>%
    rbind(hosp1$ageMid)
hosp_preds <- exp(samples %*% hosp_preds)
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
    mutate(ageMid = hosp1$ageMid)
hosp_preds1_sum <- hosp_preds1 %>%
    apply(2, function(x) {
        tibble(mean = mean(x), LCI = quantile(x, probs = 0.005), UCI = quantile(x, probs = 0.995))
    }) %>%
    bind_rows() %>%
    mutate(ageMid = hosp1$ageMid)

## fitted plot
p <- ggplot(hosp, aes(x = ageMid)) +
    geom_point(aes(y = prop)) +
    geom_line(aes(y = mean), data = hosp_preds_sum) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), data = hosp_preds_sum, alpha = 0.5) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), data = hosp_preds1_sum, alpha = 0.5) +
    xlab("Age") + ylab("Probability of hospitalisation")
ggsave("fittedHosp.svg", p)

## pick closest runs to data
closeRuns <- apply(hosp_preds, 1, function(x, p) {
        sum(abs(x - p))
    }, p = hosp1$prop) %>%
    sort.list() %>%
    {hosp_preds[.[1:10], ]} %>%
    as_tibble() %>%
    mutate(run = as.character(1:n())) %>%
    gather(ageMid, prop, -run) %>%
    mutate(ageMid = gsub("V", "", ageMid)) %>%
    mutate(ageMid = hosp1$ageMid[as.numeric(ageMid)])

## check that some are similar to data
ggplot(hosp, aes(x = ageMid)) +
    geom_point(aes(y = prop)) +
    geom_line(aes(y = mean), data = hosp_preds_sum) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), data = hosp_preds_sum, alpha = 0.5) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), data = hosp_preds1_sum, alpha = 0.5) +
    geom_line(aes(y = prop, colour = run), data = closeRuns) +
    guides(colour = FALSE) +
    xlab("Age") + ylab("Probability of hospitalisation")

## posterior correlation
cor(samples)

## subsample posterior
samples <- samples[sample.int(nrow(samples), 30000), ]

## fit range of finite mixture models
mod <- densityMclust(samples)

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

## plot posterior and approximation
p <- as_tibble(samples) %>%
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
ggsave("mixturePosterior.png", p)

## save mclust object
saveRDS(mod, "hospitalisation.rds")
