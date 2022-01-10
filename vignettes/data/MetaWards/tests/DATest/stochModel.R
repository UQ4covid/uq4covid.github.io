## load libraries
library(tidyverse)
library(Rcpp)

## create output directory
dir.create("outputs")

## read in NGM function
source("NGM.R")

## read in parameters
pars <- read_delim("disease.dat", delim = " ") %>%
    select(ends_with("_8"), nu = `beta[1]`, nuA = `beta[6]` / `beta[1]`, output)
inputs <- readRDS("inputs.rds")
colnames(pars) <- gsub("_8", "", colnames(pars))
colnames(pars) <- gsub("\\.", "", colnames(pars))
pars <- pars %>%
    select(nuold = nu, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD, output) %>%
    inner_join(
        mutate(inputs, gammaE = 1 / TE, gammaA = 1 / (TP + TI1 + TI2), gammaP = 1 / TP, gammaI1 = 1 / TI1, gammaI2 = 1 / TI2), 
    by = "output") %>%
    mutate(nu = pmap_dbl(list(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2), 
        function(R0, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2, S0, N) {
            NGM(R0 = R0, nu = NA, S0, N, nuA, gammaE, pEP, gammaA, gammaP, gammaI1, pI1H, pI1D, gammaI2)$nu
        }, S0 = 149, N = 150))

## plot old scaling against new
p <- ggplot(pars) +
    geom_point(aes(x = nuold, y = nu)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed")
ggsave("outputs/nuscaling.pdf", p)

## extract parameters for simulation   
pars <- select(slice(pars, 6), nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD) %>%
    unlist()

## set initial counts
u <- numeric(12)
u[1] <- 149
u[2] <- 1

## set seed
set.seed(4578)

## try discrete-time model
sourceCpp("discreteStochModel.cpp")
disSims <- list()
for(i in 1:50) {
    disSims[[i]] <- discreteStochModel(pars, 0, 150, u)
}
stageNms <- c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH")
disSims <- map(disSims, ~as_tibble(.)) %>%
    bind_rows(.id = "rep") %>%
    set_names(c("rep", "t", stageNms)) %>%
    mutate(t = t + 1)

## extract simulation closest to median
medRep <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        median = median(n),
        .groups = "drop"
    ) %>%
    inner_join(
        pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n"),
        by = c("t", "var")
    ) %>%
    mutate(diff = abs(median - n)) %>%
    group_by(rep) %>%
    summarise(diff = sum(diff), .groups = "drop") %>%
    arrange(diff) %>%
    slice(2) %>%
    inner_join(disSims, by = "rep") %>%
    select(!c(diff, rep))

## apply underdispersion model to counts
obsScale <- 0.8
scaleFn <- function(count, obsScale) {
    ## create incidence over time
    inc <- diff(c(0, count))
    inc <- rep(1:length(count), times = inc)
    inc <- sample(inc, round(length(inc) * obsScale))
    inc <- table(inc)
    count <- numeric(length(count))
    count[as.numeric(names(inc))] <- inc
    cumsum(count)
}
medRep <- mutate(medRep, DIobs = scaleFn(DI, obsScale = obsScale)) %>%
    mutate(medRep, DHobs = scaleFn(DH, obsScale))

## plot replicates
p <- pivot_longer(disSims, !c(rep, t), names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    ggplot(aes(x = t)) +
        geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
        geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
        geom_line(aes(y = median)) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(select(medRep, -ends_with("obs")), !t, names_to = "var", values_to = "n"),
            col = "red", linetype = "dashed"
        ) +
        geom_line(
            aes(y = n), 
            data = pivot_longer(select(medRep, t, ends_with("obs")), !t, names_to = "var", values_to = "n") %>%
                mutate(var = gsub("obs", "", var)),
            col = "blue", linetype = "dashed"
        ) +
        facet_wrap(~var, scales = "free") +
        xlab("Days") + 
        ylab("Counts")
ggsave("outputs/sims.pdf", p)

saveRDS(medRep, "outputs/disSims.rds")
saveRDS(pars, "outputs/pars.rds")

