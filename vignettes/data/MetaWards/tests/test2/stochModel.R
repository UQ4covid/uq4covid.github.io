## load libraries
library(tidyverse)
library(Rcpp)
library(SimBIID)

## set up simulation model
transitions <- c(
    "S -> beta * S * (P + I1) / (S + E + P + I1 + DI) -> E", 
    "E -> gammaE * E -> P", 
    "P -> gammaP * P -> I1", 
    "I1 -> gammaI1 * I1 -> DI"
)
compartments <- c("S", "E", "P", "I1", "DI")
parnms <- c("beta", "gammaE", "gammaP", "gammaI1")
model <- mparseRcpp(
    transitions = transitions, 
    compartments = compartments,
    pars = parnms,
    tspan = TRUE
)

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
names(pars) <- parnms

## set population size
N <- round(pluck(read_csv("age_seeds.csv", col_names = FALSE), "X2") * 100000)

## compile and run model
sims <- run(
    model = model,
    pars = pars,
    tstart = 1,
    tstop = 200,
    tspan = 2:200,
    nrep = 50,
    u = c(S = N[1] - 10, E = 10, P = 0, I1 = 0, DI = 0)
)

## try discrete-time model
sourceCpp("discreteStochModel.cpp")
disSims <- list()
for(i in 1:10) {
    disSims[[i]] <- discreteStochModel(pars, 1, 200, c(S = N[1] - 10, E = 10, P = 0, I1 = 0, DI = 0))
}

