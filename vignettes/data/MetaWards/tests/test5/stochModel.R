## load libraries
library(tidyverse)
library(Rcpp)

## read in parameters
pars <- read_delim("disease.dat", delim = " ") %>%
    select(ends_with("_1"), nu = `beta[1]`, nuA = `beta[6]`)
colnames(pars) <- gsub("_1", "", colnames(pars))
colnames(pars) <- gsub("\\.", "", colnames(pars))
pars <- pars %>%
    select(nu, nuA, pE, pEP, pA, pP, pI1, pI1H, pI1D, pI2, pH, pHD) %>%
    unlist()

## set initial counts
u <- matrix(0, 12, 8)
u[1, ] <- N
u[1, 1] <- N[1] - 10
u[2, 1] <- 10

## set contact matrices
contact1 <- read_csv("contact_matrix1.csv", col_names = FALSE) %>%
    as.matrix()
contact2 <- read_csv("contact_matrix2.csv", col_names = FALSE) %>%
    as.matrix()
    
## set up lockdown time
locktime <- as.Date("23/03/2020", format = "%d/%m/%Y")
locktime <- locktime - as.Date("01/01/2020", format = "%d/%m/%Y")
locktime <- as.numeric(locktime) - 1
print(paste0("locktime = ", locktime))

## try discrete-time model
sourceCpp("discreteStochModel.cpp")
disSims <- list()
for(i in 1:50) {
    disSims[[i]] <- discreteStochModel(pars, 150, u, contact1, contact2, locktime)
}
stageNms <- map(c("S", "E", "A", "RA", "P", "I1", "DI", "I2", "RI", "H", "RH", "DH"), ~paste0(., 1:8)) %>%
    reduce(c)
disSims <- map(disSims, ~as_tibble(.)) %>%
    bind_rows(.id = "rep") %>%
    set_names(c("rep", "t", stageNms)) %>%
    mutate(t = t + 1)

