## load libraries
library(mclust)
library(tidyverse)

## read in FMM
model <- readRDS("hospitalisation.rds")

## produce large number of samples from model
nsamps <- 1000
sims <- sim(model$modelName, model$parameters, nsamps)[, -1]

## set range of intercepts
alphaEP <- seq(-10, -2.5, length.out = 9)
alphaID <- seq(-10, -2.5, length.out = 9)
alphaHD <- seq(-10, -2.5, length.out = 9)

## expand to grid of values
alphas <- expand.grid(alphaEP, alphaID, alphaHD) %>%
    set_names(c("alphaEP", "alphaID", "alphaHD"))

## set number of individuals in each age-class
nages <- rep(10000, 8)

## read in age ranges
ageMid <- c(4.5, 14.5, 24.5, 34.5, 44.5, 54.5, 64.5, 75.0)

## set up slopes
ageSlopes <- sims[, 2] %*% t(ageMid)

## generate probs for each combination of intercepts
probs <- apply(alphas, 1, function(a, sl) {
    list(
        EP = exp(a[1] + sl),
        ID = exp(a[2] + sl),
        HD = exp(a[3] + sl)
    )
}, sl = ageSlopes)

## extract probabilities for different transitions
probsEP <- map(probs, 1)
probsID <- map(probs, 2)
probsHD <- map(probs, 3)

## generate IH probabilities from posterior draws
probsIH <- t(apply(sims, 1, function(x, ageMid) {
    exp(x[1] + x[2] * ageMid)
}, ageMid = ageMid))

## now do movements through pathways
## (rmultinom doesn't vectorise, so this takes a while)

## pmap loops over combinations
nHD <- pmap(list(probsEP, probsID, probsHD), function(EP, ID, HD, IH, nages) {
    nHD <- matrix(NA, nrow(EP), length(nages))
    ## loop over ages
    for(i in 1:ncol(EP)) {
        nP <- rbinom(rep(1, nrow(EP)), rep(nages[i], nrow(EP)), EP[, i])
        ## loop over MC samples
        temp <- pmap(list(ID[, i], IH[, i], nP), function(x, y, n) {
            if((x + y) > 1 | is.na(n)) {
                return(c(NA, NA, NA))
            } else {
                return(rmultinom(1, n, c(x, y, 1 - (x + y))))
            }
        })
        nID <- map_int(temp, 1)
        nIH <- map_int(temp, 2)
        nHD[!is.na(nIH), i] <- rbinom(rep(1, sum(!is.na(nIH))), nIH[!is.na(nIH)], HD[!is.na(nIH), i])
    }
    nHD / nages[1]
}, IH = probsIH, nages = nages)

## summarise progression and bind to inputs
nHDsum <- map(nHD, ~apply(., 2, function(x) {
        c(mean(x, na.rm = TRUE), 
          quantile(x, probs = 0.025, na.rm = TRUE), 
          quantile(x, probs = 0.975, na.rm = TRUE), 
          sum(is.na(x)))
    })) %>%
    map(~{
        t(.) %>%
        as.data.frame() %>%
        set_names(c("mean", "LCI", "UCI", "nmiss")) %>%
        mutate(age = 1:nrow(.))
    }) %>%
    bind_rows(.id = "input") %>%
    mutate(input = as.numeric(input)) %>%
    mutate(age = ageMid[age]) %>%
    left_join(mutate(alphas, input = 1:n()), by = "input") %>%
    mutate(pmiss = nmiss / nsamps)

## plot summaries
ggplot(nHDsum, aes(x = age, y = mean, group = input, colour = pmiss)) +
    geom_line() +
    facet_wrap(~alphaEP)
ggplot(nHDsum, aes(x = age, y = mean, group = input, colour = pmiss)) +
    geom_line() +
    facet_wrap(~alphaID)
ggplot(nHDsum, aes(x = age, y = mean, group = input, colour = pmiss)) +
    geom_line() +
    facet_wrap(~alphaHD)

## extract older age class
temp <- map(nHD, ~tibble(nHD = .[, 8]))
temp1 <- as_tibble(alphas)
temp1$nHD <- temp

mutate(temp1, pmiss = map_dbl(nHD, ~sum(is.na(.$nHD))) / nsamps) %>%
    select(-nHD) %>%
    group_by(alphaEP, alphaID) %>%
    summarise(pmiss = mean(pmiss), .groups = "drop") %>%
    ggplot(aes(x = alphaEP, y = alphaID, fill = pmiss)) +
    geom_tile()

mutate(temp1, p = map_dbl(nHD, ~mean(.$nHD, na.rm = TRUE))) %>%
    select(-nHD) %>%
    group_by(alphaEP, alphaID) %>%
    summarise(p = max(p), .groups = "drop") %>%
    ggplot(aes(x = alphaEP, y = alphaID, fill = p)) +
    geom_tile() +
    scale_fill_viridis_c()
