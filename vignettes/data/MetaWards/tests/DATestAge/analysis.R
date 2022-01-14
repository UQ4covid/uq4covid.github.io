## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)

## source simulation model
sourceCpp("discreteStochModel.cpp")

## read in truncated Skellam sampler
source("trSkellam.R")

## source function to run PF and return log-likelihood
source("PF.R")

## read in simulated data and generate incidence curves
data <- readRDS("outputs/disSims.rds")

## read in parameters, remove guff and reorder
pars <- read_delim("disease.dat", delim = " ") %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta"), repeats, starts_with(".lock"), .p_home_weekend)) %>%
    select(nu, nuA, !output)

## read in contact matrix
contact <- read_csv("POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## solution to round numbers preserving sum
## adapted from:
## https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
    y <- floor(x)
    indices <- tail(order(x - y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
}

## set up number of initial individuals in each age-class
N <- 150
N <- smart_round(read_csv("age_seeds.csv", col_names = FALSE)$X2 * N)
I0 <- smart_round(read_csv("age_seeds.csv", col_names = FALSE)$X2 * 1)
N <- N - I0

## set initial counts
u <- matrix(0, 12, 8)
u[1, ] <- N
u[2, ] <- I0

## set seed for reproducibility
set.seed(666)

## run model with no model discrepancy
runs_nomd <- PF(pars, C = contact, data = data, u = u, ndays = 30, npart = 100, obsScale = 0.8, MD = FALSE)

## repeat but adding some model discrepancy
runs_md <- PF(pars, C = contact, data = data, u = u, ndays = 30, npart = 100, obsScale = 0.8, MD = TRUE, disScale = 0.1)

## run model with no model discrepancy and throw out sims corresponding to true parameters
## if you want to save out some runs, you can only run for a single design point at a time
## due to parallelisation - I could fix, but not right now
runs_nomd <- PF(pars[6, ], C = contact, data = data, u = u, ndays = 30, npart = 100, obsScale = 0.8, MD = FALSE, whichSave = 1)

## plot particle estimates of states (unweighted)
sims_nomd <- map(sims, ~map(., ~as.vector(t(.)))) %>%
    map(~do.call("rbind", .)) %>%
    map(as_tibble) %>%
    bind_rows(.id = "t") %>%
    mutate(t = as.numeric(t))
stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., 1:8)) %>%
    reduce(c)
colnames(sims_nomd) <- c("t", stageNms)

## plot replicates
p <- pivot_longer(sims_nomd, !t, names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(age = gsub("[^0-9]", "", var)) %>%
    mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    ggplot(aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
    geom_line(
        aes(y = n),
        data = pivot_longer(select(data, !ends_with("obs")) %>% filter(t <= 30), !t, names_to = "var", values_to = "n") %>%
            mutate(age = gsub("[^0-9]", "", var)) %>%
            mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
            mutate(var = gsub("one", "1", var)) %>%
            mutate(var = gsub("two", "2", var)),
        col = "red", linetype = "dashed"
    ) +
    geom_line(
        aes(y = n),
        data = pivot_longer(select(data, t, ends_with("obs")) %>% filter(t <= 30), !t, names_to = "var", values_to = "n") %>%
            mutate(var = gsub("obs", "", var)) %>%
            mutate(age = gsub("[^0-9]", "", var)) %>%
            mutate(var = gsub("[^a-zA-Z]", "", var)),
        col = "blue", linetype = "dashed"
    ) +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts")
ggsave("sims_nomd.pdf", p, width = 10, height = 10)

## repeat but adding some model discrepancy
runs_md <- PF(pars[6, ], C = contact, data = data, u = u, ndays = 30, npart = 100, obsScale = 0.8, MD = TRUE, disScale = 0.1, whichSave = 1)

## plot particle estimates of states (unweighted)
sims_md <- map(sims, ~map(., ~as.vector(t(.)))) %>%
    map(~do.call("rbind", .)) %>%
    map(as_tibble) %>%
    bind_rows(.id = "t") %>%
    mutate(t = as.numeric(t))
stageNms <- map(c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"), ~paste0(., 1:8)) %>%
    reduce(c)
colnames(sims_md) <- c("t", stageNms)

## plot replicates
p <- pivot_longer(sims_md, !t, names_to = "var", values_to = "n") %>%
    group_by(t, var) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(age = gsub("[^0-9]", "", var)) %>%
    mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    ggplot(aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ), alpha = 0.5) +
    geom_line(aes(y = median)) +
    geom_line(
        aes(y = n),
        data = pivot_longer(select(data, !ends_with("obs")) %>% filter(t <= 30), !t, names_to = "var", values_to = "n") %>%
            mutate(age = gsub("[^0-9]", "", var)) %>%
            mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
            mutate(var = gsub("one", "1", var)) %>%
            mutate(var = gsub("two", "2", var)),
        col = "red", linetype = "dashed"
    ) +
    geom_line(
        aes(y = n),
        data = pivot_longer(select(data, t, ends_with("obs")) %>% filter(t <= 30), !t, names_to = "var", values_to = "n") %>%
            mutate(var = gsub("obs", "", var)) %>%
            mutate(age = gsub("[^0-9]", "", var)) %>%
            mutate(var = gsub("[^a-zA-Z]", "", var)),
        col = "blue", linetype = "dashed"
    ) +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts")
ggsave("sims_md.pdf", p, width = 10, height = 10)

## plot both together
p <- mutate(sims_nomd, type = "No MD") %>%
    rbind(mutate(sims_md, type = "MD")) %>%
    pivot_longer(!c(t, type), names_to = "var", values_to = "n") %>%
    group_by(t, var, type) %>%
    summarise(
        LCI = quantile(n, probs = 0.025),
        LQ = quantile(n, probs = 0.25),
        median = median(n),
        UQ = quantile(n, probs = 0.75),
        UCI = quantile(n, probs = 0.975),
        .groups = "drop"
    ) %>%
    mutate(age = gsub("[^0-9]", "", var)) %>%
    mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
    mutate(var = gsub("one", "1", var)) %>%
    mutate(var = gsub("two", "2", var)) %>%
    ggplot(aes(x = t)) +
    geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = type), alpha = 0.5) +
    geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = type), alpha = 0.5) +
    geom_line(aes(y = median, colour = type)) +
    geom_line(
        aes(y = n),
        data = pivot_longer(select(data, !ends_with("obs")) %>% filter(t <= 30), !t, names_to = "var", values_to = "n") %>%
            mutate(age = gsub("[^0-9]", "", var)) %>%
            mutate(var = gsub("[^a-zA-Z]", "", var)) %>%
            mutate(var = gsub("one", "1", var)) %>%
            mutate(var = gsub("two", "2", var)),
        col = "red", linetype = "dashed"
    ) +
    geom_line(
        aes(y = n),
        data = pivot_longer(select(data, t, ends_with("obs")) %>% filter(t <= 30), !t, names_to = "var", values_to = "n") %>%
            mutate(var = gsub("obs", "", var)) %>%
            mutate(age = gsub("[^0-9]", "", var)) %>%
            mutate(var = gsub("[^a-zA-Z]", "", var)),
        col = "blue", linetype = "dashed"
    ) +
    facet_grid(var ~ age, scales = "free") +
    xlab("Days") + 
    ylab("Counts")
ggsave("sims_combined.pdf", p, width = 10, height = 10)
