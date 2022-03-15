## load libraries
library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(sitmo)

## source Rcpp PF code
sourceCpp("PF.cpp")

## source function to run PF and return log-likelihood
source("PF.R")

## set wave
wave <- 1

## read in simulated data and generate incidence curves
data <- readRDS("outputs/disSims.rds")

## read in parameters, remove guff and reorder
pars <- readRDS(paste0("wave", wave, "/disease.rds")) %>%
    rename(nu = `beta[1]`, nuA = `beta[6]`) %>%
    select(!c(starts_with("beta"), repeats)) %>%
    select(nu, nuA, !output)

## read in contact matrix
contact <- read_csv("inputs/POLYMOD_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## read in initial conditions
u1 <- readRDS("outputs/u1.rds")
u1_moves <- readRDS("outputs/u1_moves.rds")

## run PF with some model discrepancy
runs_md <- runs_md <- PF(pars[1:8, ], C = contact, data = data, u1_moves = u1_moves,
    u1 = u1, ndays = 50, npart = 10, MD = TRUE, a_dis = 0.05, b_dis = 0.05, saveAll = NA)

## save outputs
saveRDS(runs_md, paste0("wave", wave, "/runs_md.rds"))
