## load libraries
library(mclust)
library(tidyverse)
library(nimble)
library(coda)
library(GGally)
library(svglite)

## load data
hosp <- read_csv("hospitalisationImperial.csv", col_names = TRUE) %>%
    rename(ageCat = X1, severe = `Severe cases`, all = `All cases`) %>%
    select(ageCat, severe, all) %>%
    mutate(ageCat = gsub("–", "-", ageCat)) %>%
    mutate(ageCat = gsub(" years", "", ageCat))

## match data roughly to age-classes used in the model
hosp <- mutate(hosp, ageCat = ifelse(ageCat == "70-79", "70+", ageCat)) %>%
    mutate(ageCat = ifelse(ageCat == "≥80", "70+", ageCat)) %>%
    group_by(ageCat) %>%
    summarise_all(sum) %>%
    arrange(ageCat)

## extract mid-points of age-classes
hosp <- hosp %>%
    mutate(ageCat = ifelse(ageCat == "70+", "70-80", ageCat)) %>%
    separate(ageCat, c("LB", "UB"), sep = "-", remove = FALSE) %>%
    mutate_at(vars(c("UB", "LB")), as.numeric) %>%
    mutate(ageMid = (LB + UB) / 2) %>%
    dplyr::select(-LB, -UB) %>%
    mutate(pS = severe / all)

## fit Bayesian model
code <- nimbleCode({
    
    ## fit log-linear model to 1 / rates
    for(k in 1:nAge) {
        severe[k] ~ dbin(pS[k], all[k])
        log(pS[k]) <- alpha + eta * ageMid[k]
    }
    
    ## priors
    alpha ~ dnorm(0, 5)
    eta ~ T(dnorm(0, 1), 0, )
})

## set up other components of model
consts <- list(
    nAge = nrow(hosp),
    ageMid = hosp$ageMid,
    all = hosp$all)
data <- list(severe = hosp$severe)

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
pdf("traces.pdf")
plot(run$samples)
dev.off()
samples <- as.matrix(run$samples)

## predictions
hosp_preds <- matrix(1, nrow = 1, ncol = nrow(hosp)) %>%
    rbind(hosp$ageMid)
hosp_preds <- exp(samples %*% hosp_preds)
hosp_preds1 <- rbind(hosp$all, hosp_preds) %>%
    apply(2, function(x) {
        rbinom(rep(1, length(x) - 1), x[1], x[-1]) / x[1]
    })
hosp_preds <- apply(hosp_preds, 2, function(x) {
    tibble(mean = mean(x), LCI = quantile(x, probs = 0.005), UCI = quantile(x, probs = 0.995))
}) %>%
    bind_rows() %>%
    mutate(ageMid = hosp$ageMid)
hosp_preds1 <- apply(hosp_preds1, 2, function(x) {
    tibble(mean = mean(x), LCI = quantile(x, probs = 0.005), UCI = quantile(x, probs = 0.995))
}) %>%
    bind_rows() %>%
    mutate(ageMid = hosp$ageMid)

## fitted plot
p <- ggplot(hosp, aes(x = ageMid)) +
    geom_point(aes(y = pS)) +
    geom_line(aes(y = mean), data = hosp_preds) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), data = hosp_preds, alpha = 0.5) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), data = hosp_preds1, alpha = 0.5) +
    xlab("Age") + ylab("Probability of hospitalisation")
ggsave("fittedLnBayes.pdf", p)
ggsave("fittedLnBayes.svg", p)

## posterior correlation
cor(samples)

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
ggsave("mixturePosterior.pdf", p)
ggsave("mixturePosterior.svg", p)
ggsave("mixturePosterior.png", p)

## save mclust object
saveRDS(mod, "hospitalisation.rds")
