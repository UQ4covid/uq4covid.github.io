## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)

## source simulation model
sourceCpp("discreteStochModel.cpp")

## read in NGM function
source("NGM.R")

## read in truncated Skellam sampler
source("trSkellam.R")

## source function to run PF and return log-likelihood
source("PF.R")

## read in simulated data and generate incidence curves
data <- readRDS("outputs/disSims.rds")

## read in parameters and rescale nu for a single population
pars <- read_delim("disease.dat", delim = " ") %>%
    select(ends_with("_8"), nu = `beta[1]`, nuA = `beta[6]` / `beta[1]`, output)
inputs <- readRDS("inputs.rds")
colnames(pars) <- gsub("_8", "", colnames(pars))
colnames(pars) <- gsub("\\.", "", colnames(pars))
pars <- pars %>%
    select(pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD, output) %>%
    inner_join(
        mutate(inputs, gammaE = 1 / TE, gammaA = 1 / (TP + TI1 + TI2), gammaP = 1 / TP, gammaI1 = 1 / TI1, gammaI2 = 1 / TI2), 
    by = "output") %>%
    mutate(nu = pmap_dbl(list(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2), 
        function(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2, S0, N) {
            NGM(R0 = R0, nu = NA, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2)$nu
        }, S0 = 149, N = 150))
        
## select parameter sets
pars <- select(pars, nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD)

## run particle filter over 1st 10 days and calculate log-likelihood
u <- numeric(12)
u[1] <- 149
u[2] <- 1

## set seed for reproducibility
set.seed(666)

## run model with no model discrepancy
runs_nomd <- PF(pars, data = data, u = u, ndays = 30, npart = 100, obsScale = 0.8, MD = FALSE)

## repeat but adding some model discrepancy
runs_md <- PF(pars, data = data, u = u, ndays = 30, npart = 100, obsScale = 0.8, MD = TRUE, disScale = 0.1)

## run model with no model discrepancy and throw out sims corresponding to true parameters
## if you want to save out some runs, you can only run for a single design point at a time
## due to parallelisation - I could fix, but not right now
runs_nomd <- PF(pars[6, ], data = data, u = u, ndays = 30, npart = 100, obsScale = 0.8, MD = FALSE, whichSave = 1)

## plot particle estimates of states (unweighted)
sims_nomd <- map(sims, as_tibble) %>%
    bind_rows(.id = "t") %>%
    mutate(t = factor(t, levels = sort(as.numeric(unique(t)))))
colnames(sims_nomd) <- c("t", "S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
p <- pivot_longer(sims_nomd, !t, names_to = "var", values_to = "count") %>%
    ggplot(aes(x = t, y = count)) +
    geom_violin() +
    geom_point(data = pivot_longer(
        select(data, !ends_with("obs")) %>% 
            filter(t <= 30), 
        !t, names_to = "var", values_to = "count")) +
    geom_point(data = pivot_longer(
        select(data, t, ends_with("obs")) %>% 
            filter(t <= 30), 
        !t, names_to = "var", values_to = "count") %>%
            mutate(var = gsub("obs", "", var)),
        col = "blue") +
    facet_wrap(~var, scales = "free")
ggsave("sims_nomd.pdf", p)

## repeat but adding some model discrepancy
runs_md <- PF(pars[6, ], data = data, u = u, ndays = 30, npart = 100, obsScale = 0.8, MD = TRUE, disScale = 0.1, whichSave = 1)

## plot particle estimates of states (unweighted)
sims_md <- map(sims, as_tibble) %>%
    bind_rows(.id = "t") %>%
    mutate(t = factor(t, levels = sort(as.numeric(unique(t)))))
colnames(sims_md) <- c("t", "S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
p <- pivot_longer(sims_md, !t, names_to = "var", values_to = "count") %>%
    ggplot(aes(x = t, y = count)) +
    geom_violin() +
    geom_point(data = pivot_longer(
            select(data, !ends_with("obs")) %>% 
            filter(t <= 30), 
        !t, names_to = "var", values_to = "count")) +
    geom_point(data = pivot_longer(
                select(data, t, ends_with("obs")) %>% 
                filter(t <= 30), 
            !t, names_to = "var", values_to = "count") %>%
        mutate(var = gsub("obs", "", var)),
        col = "blue") +
    facet_wrap(~var, scales = "free")
ggsave("sims_md.pdf", p)

## plot both together
p <- mutate(sims_nomd, type = "No MD") %>%
    rbind(mutate(sims_md, type = "MD")) %>%
    pivot_longer(!c(t, type), names_to = "var", values_to = "count") %>%
    ggplot(aes(x = t, y = count)) +
    geom_violin(aes(colour = type)) +
    geom_point(data = pivot_longer(
        select(data, !ends_with("obs")) %>% 
            filter(t <= 30), 
        !t, names_to = "var", values_to = "count")) +
    geom_point(data = pivot_longer(
        select(data, t, ends_with("obs")) %>% 
            filter(t <= 30), 
        !t, names_to = "var", values_to = "count") %>%
            mutate(var = gsub("obs", "", var)),
        col = "blue") +
    facet_wrap(~var, scales = "free")
ggsave("sims_combined.pdf", p)        
