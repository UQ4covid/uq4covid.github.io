## load libraries
library(mclust)
library(MASS)
library(tidyverse)
library(GGally)

## set seed
set.seed(666)

## load FMM objects
pathways <- readRDS("../../inputs/pathways.rds")

## generate random samples and grab threshold from these
nsamp <- 1000000
sims <- sim(pathways$modelName, pathways$parameters, nsamp)[, -1]
temp <- dens(sims, pathways$modelName, pathways$parameters, logarithm = TRUE)
sims <- as.data.frame(sims) %>%
    set_names(c("alphaEP", "alphaID", "alphaHD", "alphaIH", "eta")) %>%
    mutate(dens = temp) %>%
    filter(eta > 0 & eta < 1) %>%
    filter(alphaEP > -20 & alphaEP < 0) %>%
    filter(alphaID > -20 & alphaID < 0) %>%
    filter(alphaHD > -20 & alphaHD < 0) %>%
    filter(alphaIH > -20 & alphaIH < 0)

## generate p.d.f. threshold
paththresh <- round(quantile(sims$dens, probs = 0.001), 3)

## produce plot of valid samples
nsamp <- 1000000
sims <- data.frame(
        alphaEP = runif(nsamp, -20, 0),
        alphaID = runif(nsamp, -20, 0),
        alphaHD = runif(nsamp, -20, 0),
        alphaIH = runif(nsamp, -20, 0),
        eta = runif(nsamp, 0, 1)
    )
sims$dens <- dens(as.matrix(sims), pathways$modelName, pathways$parameters, logarithm = TRUE)
sims <- filter(sims, dens > paththresh)
while(nrow(sims) < 1000) {
    sims1 <- data.frame(
        alphaEP = runif(nsamp, -20, 0),
        alphaID = runif(nsamp, -20, 0),
        alphaHD = runif(nsamp, -20, 0),
        alphaIH = runif(nsamp, -20, 0),
        eta = runif(nsamp, 0, 1)
    )
    sims1$dens <- dens(as.matrix(sims1), pathways$modelName, pathways$parameters, logarithm = TRUE)
    sims1 <- filter(sims1, dens > paththresh)
    sims <- rbind(sims, sims1)
}

## plot valid samples
p <- ggpairs(sims, columns = 1:5, upper = "blank") 
ggsave("pathwaysLimits.pdf", p)
saveRDS(paththresh, "../../inputs/pathThresh.rds")

