## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(abind)
library(sitmo)

## read in truncated Skellam sampler
source("trSkellam.R")

## source simulation model
sourceCpp("PF.cpp")

## pre-emptive garbage collect
gc()

## source function to run PF and return log-likelihood
source("PF.R")

## read in simulated data and generate incidence curves
data <- readRDS("outputs/disSims.rds")

## read in parameters, remove guff and reorder
pars <- read_delim("disease.dat", delim = " ") %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta"), repeats, starts_with(".lock"), .p_home_weekend)) %>%
    select(nu, nuA, !output) %>%
    as.data.frame()

## read in contact matrix
contact <- read_csv("POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## read in initial conditions
u1 <- readRDS("outputs/u1.rds")
u1_moves <- readRDS("outputs/u1_moves.rds")

## set seed for reproducibility
set.seed(666)

# ## run model with no model discrepancy
# runs_nomd <- PF(pars[1:10, ], C = contact, data = data, u = u, ndays = 30, npart = 100, MD = FALSE, saveAll = NA)

# ## repeat but adding some model discrepancy
# runs_md <- PF(pars, C = contact, data = data, u = u, ndays = 30, npart = 100, MD = TRUE, a_dis = 0.05, b_dis = 0.05, saveAll = NA)

 ## set up plot data
 plot_data <- pivot_longer(filter(data, t <= 50), !t, names_to = "var", values_to = "n") %>%
     mutate(age = gsub('^(?:[^_]*_)(.*)', '\\1', var)) %>%
     mutate(LAD = gsub('^(?:[^_]*_)(.*)', '\\1', age)) %>%
     mutate(age = gsub('(.*)_[0-9]*', '\\1', age)) %>%
     mutate(obs = grepl("obs", var)) %>%
     mutate(var = gsub('^(.*)_[0-9]*_.*', '\\1', var)) %>%
     group_by(t, var, age, obs) %>%
     summarise(n = sum(n), .groups = "drop") %>%
     mutate(var = gsub("one", "1", var)) %>%
     mutate(var = gsub("two", "2", var)) %>%
     mutate(age = gsub("obs", "", age))
 
 for(k in 246) {
     ## run model with no model discrepancy and throw out sims corresponding to true parameters
     ## if you want to save out some runs, you can only run for a single design point at a time
     ## due to parallelisation - I could fix, but not right now
     runs_nomd <- PF(pars[k, ], C = contact, data = data, u1_moves = u1_moves,
                     u1 = u1, ndays = 50, npart = 10, MD = FALSE, saveAll = TRUE)
 
     ## plot particle estimates of states (unweighted)
     sims_nomd <- map(runs_nomd$particles[[1]], ~{
             map(., ~{
                 x <- apply(., c(1, 2), sum)
                 colnames(x) <- paste0("age", 1:ncol(x))
                 as_tibble(x) %>%
                     mutate(var = c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"))
             }) %>%
             bind_rows(.id = "t")
         }) %>%
         bind_rows(.id = "particle") %>%
         pivot_longer(!c(particle, t, var), names_to = "age", values_to = "n") %>%
         mutate(age = as.numeric(gsub("age", "", age))) %>%
         mutate(t = as.numeric(t) - 1) %>%
         mutate(var = gsub("one", "1", var)) %>%
         mutate(var = gsub("two", "2", var))

     ## repeat but adding some model discrepancy
     runs_md <- PF(pars[k, ], C = contact, data = data, u1_moves = u1_moves,
                   u1 = u1, ndays = 50, npart = 10, MD = TRUE, a_dis = 5e-10, b_dis = 0, saveAll = TRUE)
 
     ## plot particle estimates of states (unweighted)
     sims_md <- map(runs_md$particles[[1]], ~{
             map(., ~{
                 x <- apply(., c(1, 2), sum)
                 colnames(x) <- paste0("age", 1:ncol(x))
                 as_tibble(x) %>%
                     mutate(var = c("S", "E", "A", "RA", "P", "Ione", "DI", "Itwo", "RI", "H", "RH", "DH"))
             }) %>%
                 bind_rows(.id = "t")
         }) %>%
         bind_rows(.id = "particle") %>%
         pivot_longer(!c(particle, t, var), names_to = "age", values_to = "n") %>%
         mutate(age = as.numeric(gsub("age", "", age))) %>%
         mutate(t = as.numeric(t) - 1) %>%
         mutate(var = gsub("one", "1", var)) %>%
         mutate(var = gsub("two", "2", var))
         
     ## plot both together nationally
     p <- mutate(sims_nomd, type = "No MD") %>%
         rbind(mutate(sims_md, type = "MD")) %>%
         group_by(t, var, type) %>%
         summarise(
             LCI = quantile(n, probs = 0.025),
             LQ = quantile(n, probs = 0.25),
             median = median(n),
             UQ = quantile(n, probs = 0.75),
             UCI = quantile(n, probs = 0.975),
             .groups = "drop"
         ) %>%
         mutate(var = gsub("one", "1", var)) %>%
         mutate(var = gsub("two", "2", var)) %>%
         ggplot(aes(x = t)) +
         geom_ribbon(aes(ymin = LCI, ymax = UCI, fill = type), alpha = 0.5) +
         geom_ribbon(aes(ymin = LQ, ymax = UQ, fill = type), alpha = 0.5) +
         geom_line(aes(y = median, colour = type)) +
         geom_line(aes(y = n), data = filter(plot_data, !obs), col = "red", linetype = "dashed") +
         geom_line(aes(y = n), data = filter(plot_data, obs), col = "blue", linetype = "dashed") +
         facet_grid(var ~ age, scales = "free") +
         xlab("Days") +
         ylab("Counts")
     ggsave(paste0("sims_combined_", k, ".pdf"), p, width = 10, height = 10)
 }
