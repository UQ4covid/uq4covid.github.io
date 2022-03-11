## load libraries
library(mclust)
library(MASS)
library(tidyverse)
library(GGally)

## set seed
set.seed(666)

## load FMM objects
hospStays <- readRDS("../../inputs/hospStays.rds")

## generate random samples and grab threshold from these
nsamp <- 1000000
sims <- sim(hospStays$modelName, hospStays$parameters, nsamp)[, -1]
temp <- dens(sims, hospStays$modelName, hospStays$parameters, logarithm = TRUE)
sims <- as.data.frame(sims) %>%
    set_names(c("alpha", "eta")) %>%
    mutate(dens = temp) %>%
    filter(eta > 0)

## generate p.d.f. threshold 
hospthresh <- round(quantile(sims$dens, probs = 0.001), 3)

## create grid
sims <- expand.grid(alpha = seq(-2, 3, length.out = 500), eta = seq(0, 0.1, length.out = 500))
sims$dens <- dens(as.matrix(sims), hospStays$modelName, hospStays$parameters, logarithm = TRUE)
sims$valid <- ifelse(sims$dens > hospthresh, "1", "0")

## plot valid region
p <- ggplot(sims) +
    geom_raster(aes(x = alpha, y = eta, fill = valid))
ggsave("hospStaysLimits.pdf", p)
saveRDS(hospthresh, "../../inputs/hospThresh.rds")
