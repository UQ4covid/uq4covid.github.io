## load libraries
library(lhs)
library(mclust)
library(tidyverse)
library(GGally)

## source dataTools
source("R_tools/dataTools.R")

## set up parameter ranges for uniform ranges
parRanges <- data.frame(
    parameter = c("R0", "TE", "TP", "TI1", "TI2", 
                  "nuA", "lock_1_restrict", "lock_2_release"),
    lower = c(2, 0.1, 1.2, 2.8, 0.0001, 0, 0, 0),
    upper = c(4.5, 2, 3, 4.5, 0.5, 1, 1, 1),
    stringsAsFactors = FALSE
) 

## read in contact matrix to use for NGM
C <- as.matrix(read.csv("inputs/POLYMOD_matrix.csv", header = FALSE))

## generate LHS design
ndesign <- 100
design <- randomLHS(ndesign, nrow(parRanges))
colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## convert to input space
input <- convertDesignToInput(design, parRanges, "zero_one")

## generate space-filling design for other parameters

## load FMM objects
hospStays <- readRDS("inputs/hospStays.rds")
pathways <- readRDS("inputs/pathways.rds")

## generate design points for hospital stay lengths
hospStaysInput <- FMMmaximin(hospStays, ndesign, 10000) %>%
    as_tibble() %>%
    rename(alphaTH = x1, etaTH = x2)
pathwaysInput <- FMMmaximin(pathways, ndesign, 10000) %>%
    as_tibble() %>%
    rename(alphaEP = x1, alphaI1D = x2, alphaHD = x3, alphaI1H = x4, eta = x5)

## bind to design
input <- cbind(input, hospStaysInput, pathwaysInput)

## add unique hash identifier
## (at the moment don't use "a0" type ensembleID, because MetaWards
## parses to dates)
input$output <- ensembleIDGen(ensembleID = "Ens0", nrow(input))
input$repeats <- 1

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
disease <- convertInputToDisease(input, C, N, S0, ages)

## plot inputs
select(disease, ind) %>%
    mutate(valid = 1) %>%
    right_join(mutate(input, ind = 1:n()), by = "ind") %>%
    arrange(!is.na(valid)) %>%
    ggpairs(aes(colour = valid), columns = 3:17, upper = "blank")

## write text file for MetaWards
write.table(disease, "inputs/disease.dat", row.names = FALSE, sep = " ", quote = FALSE)

## save inputs for data for post-simulation runs
saveRDS(inputs, "inputs/inputs.rds")
saveRDS(parRanges, "inputs/parRanges.rds")
