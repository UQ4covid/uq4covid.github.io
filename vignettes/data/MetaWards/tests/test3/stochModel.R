## load libraries
library(tidyverse)
library(Rcpp)

## read in parameters
pars <- read_delim("disease.dat", delim = " ") %>%
    select(ends_with("_1"), nu = `beta[1]`, nuA = `beta[6]`)
colnames(pars) <- gsub("_1", "", colnames(pars))
colnames(pars) <- gsub("\\.", "", colnames(pars))
pars <- pars %>%
    select(nu, probE = pE, probP = pP, probI1 = pI1) %>%
    unlist()

## set population size
N <- c(6000, 15400, 15400, 13400, 12800, 13400, 10500, 13100)

## set initial counts
u <- matrix(0, 5, 8)
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
    disSims[[i]] <- discreteStochModel(pars, 150, u, contact)
}
stageNms <- map(c("S", "E", "P", "I1", "DI"), ~paste0(., 1:8)) %>%
    reduce(c)
disSims <- map(disSims, ~as_tibble(.)) %>%
    bind_rows(.id = "rep") %>%
    set_names(c("rep", "t", stageNms)) %>%
    mutate(t = t + 1)

