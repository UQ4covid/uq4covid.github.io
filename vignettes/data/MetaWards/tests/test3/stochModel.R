## load libraries
library(tidyverse)
library(Rcpp)

## read in parameters
pars <- read_delim("disease.dat", delim = " ") %>%
    select(ends_with("_1"), nu = `beta[1]`, nuA = `beta[6]`)
colnames(pars) <- gsub("_1", "", colnames(pars))
colnames(pars) <- gsub("\\.", "", colnames(pars))
pars <- pars %>%
    mutate_all(~ifelse(. == 1, 0.99, .)) %>%
    mutate(gammaE = -log(1 - pE)) %>%
    mutate(gammaP = -log(1 - pP)) %>%
    mutate(gammaA = -log(1 - pA)) %>%
    mutate(gammaI1 = -log(1 - pI1)) %>%
    mutate(gammaI2 = -log(1 - pI2)) %>%
    mutate(gammaH = -log(1 - pH)) %>%
    select(-pE, -pP, -pA, -pI1, -pI2, -pH) %>%
    unlist()
pars <- c(
    beta = pars["nu"], 
    gammaE = pars["gammaE"], 
    gammaP = pars["gammaP"], 
    gammaI1 = pars["gammaI1"]
)
names(pars) <- c("nu", "gammaE", "gammaP", "gammaI1")

## set population size
N <- round(pluck(read_csv("../../inputs/age_seeds.csv", col_names = FALSE), "X2") * 100000)

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
disSims <- discreteStochModel(pars, 250, u, contact)
#disSims <- list()
#for(i in 1:10) {
#    disSims[[i]] <- discreteStochModel(pars, 250, u, contact)
#}

