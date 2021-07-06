## load libraries
library(lhs)
library(mclust)
library(tidyverse)

## source dataTools
source("R_tools/dataTools.R")

## set up parameter ranges for uniform ranges
parRanges <- data.frame(
    parameter = c("R0", "TE", "TP", "TI1", "TI2", 
                  "nuA", "lock_1_restrict", "lock_2_release", "ns", "p_home_weekend"),
    lower = c(2, 0.1, 1.2, 2.8, 0.0001, 0, 0, 0, 10, 0),
    upper = c(4.5, 2, 3, 4.5, 0.5, 1, 1, 1, 100, 1),
    stringsAsFactors = FALSE
) 

## read in contact matrix to use for NGM
C <- as.matrix(read.csv("inputs/POLYMOD_matrix.csv", header = FALSE))

## generate LHS design
ndesign <- 10
design <- randomLHS(ndesign, nrow(parRanges))
colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## convert to input space
inputs <- convertDesignToInput(design, parRanges, "zero_one")

## generate space-filling design for other parameters

## load FMM objects
hospStays <- readRDS("inputs/hospStays.rds")
pathways <- readRDS("inputs/pathways.rds")

## generate design points for hospital stay lengths
hospStaysInput <- FMMmaximin(
        hospStays, 
        ndesign, 
        matrix(c(-Inf, Inf, 0, Inf), ncol = 2, byrow = TRUE)
    ) %>%
    as_tibble() %>%
    rename(alphaTH = x1, etaTH = x2)

## generate design points for other transition probabilities

## this function checks validity of inputs
pathwaysLimitFn <- function(x, ages) {
    singleProbs <- apply(x, 1, function(x, ages) {
        eta <- x[5]
        alphas <- x[-5]
        y <- sapply(alphas, function(a, eta, ages) {
            y <- exp(a + eta * ages)
            all(y >= 0 & y <= 1)
        }, eta = eta, ages = ages)
        all(y)
    }, ages = ages)
    multiProbs <- apply(x[, -c(1, 3)], 1, function(x, ages) {
        alphaI1D <- x[1]
        alphaI1H <- x[2]
        eta <- x[3]
        pI1D <- exp(alphaI1D + eta * ages)
        pI1H <- exp(alphaI1H + eta * ages)
        p <- pI1D + pI1H
        all(p >= 0 & p <= 1)
    }, ages = ages)
    multiProbs & singleProbs
}

## produces design points subject to constraints
pathwaysInput <- FMMmaximin(
        pathways, 
        ndesign,
        matrix(c(rep(c(-20, 0), times = 4), 0, 1), ncol = 2, byrow = TRUE),
        pathwaysLimitFn,
        ages = c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5)
    ) %>%
    as_tibble() %>%
    rename(alphaEP = x1, alphaI1D = x2, alphaHD = x3, alphaI1H = x4, eta = x5)

## bind to design
inputs <- cbind(inputs, hospStaysInput, pathwaysInput)

## add unique hash identifier
## (at the moment don't use "a0" type ensembleID, because MetaWards
## parses to dates)
inputs$output <- ensembleIDGen(ensembleID = "Ens0", nrow(inputs))
inputs$repeats <- 1

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
N <- smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * 56082077)
S0 <- N - smart_round(read_csv("inputs/age_seeds.csv", col_names = FALSE)$X2 * 100)
ages <- c(2.5, 11, 23.5, 34.5, 44.5, 55.5, 65.5, 75.5)

## convert input to disease
disease <- convertInputToDisease(inputs, C, N, S0, ages)
stopifnot(nrow(disease) == ndesign)

# ## plot inputs
# library(GGally)
# p <- select(inputs, -output, -repeats) %>%
#     ggpairs(upper = "blank")
# ggsave("design.pdf", p, width = 10, height = 10)

## reorder samples
inputs <- arrange(inputs, output)
disease <- arrange(disease, output)

## write text file for MetaWards
write.table(disease, "inputs/disease.dat", row.names = FALSE, sep = " ", quote = FALSE)

## save inputs for data for post-simulation runs
saveRDS(inputs, "inputs/inputs.rds")
saveRDS(parRanges, "inputs/parRanges.rds")
