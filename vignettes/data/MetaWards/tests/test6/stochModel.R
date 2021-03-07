## load libraries
library(tidyverse)
library(Rcpp)

## read in pars
pars <- read_delim("disease.dat", delim = " ") %>%
    select(-contains("beta"), -contains("lock_"), 
    nu_1 = `beta[1]`, nuA_1 = `beta[6]`, -repeats, -output) %>%
    gather(par, value) %>%
    separate(par, c("par", "age"), sep = "_") %>%
    complete(age, par) %>%
    mutate(par = gsub("\\.", "", par)) %>%
    arrange(par, age) %>%
    select(-age) %>%
    group_by(par) %>%
    fill(value) %>%
    nest() %>%
    mutate(data = map(data, "value")) %>%
    spread(par, data)

## set initial counts
u <- matrix(0, 12, 8)
u[1, ] <- N
u[1, 1] <- N[1] - 10
u[2, 1] <- 10

## set contact matrix
contact <- read_csv("contact_matrix.csv", col_names = FALSE) %>%
    as.matrix()

## try discrete-time model
sourceCpp("discreteStochModel.cpp")
disSims <- list()
for(i in 1:50) {
    disSims[[i]] <- discreteStochModel(
        pars$nu[[1]][1],
        pars$nuA[[1]][1],
        pars$pE[[1]],
        pars$pEP[[1]],
        pars$pA[[1]],
        pars$pP[[1]],
        pars$pI1[[1]],
        pars$pI1H[[1]],
        pars$pI1D[[1]],
        pars$pI2[[1]],
        pars$pH[[1]],
        pars$pHD[[1]],
        200, u, contact)
}
stageNms <- map(c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH"), ~paste0(., 1:8)) %>%
    reduce(c)
disSims <- map(disSims, ~as_tibble(.)) %>%
    bind_rows(.id = "rep") %>%
    set_names(c("rep", "t", stageNms)) %>%
    mutate(t = t + 1)

