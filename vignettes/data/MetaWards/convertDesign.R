## load libraries
library(lhs)
library(dplyr)
library(tidyr)

## source dataTools
source("R_tools/dataTools.R")

## set up parameter ranges
parRanges <- data.frame(
    parameter = c("r_zero", "latent_time", "infectious_time", "hospital_time",
                  "lock_1_restrict", "lock_2_release",
                  "alphaEA", "etaEA", 
                  "alphaIH", "etaIH", 
                  "alphaIR", "etaIR",
                  "alphaHR", "etaHR",
                  "GP_A", "GP_H"),
    lower = c(2, 4, 2, 4, 0, 0, rep(c("?", 0), 4), rep(0, 2)),
    upper = c(4, 6, 4, 12, 1, 1, rep(c("?", "?"), 4), rep(1, 2)),
    stringsAsFactors = FALSE
)

## read in contact matrix
C <- as.matrix(read.csv("inputs/contact_matrix.csv", header = FALSE))

## generate LHS design
design <- randomLHS(5, nrow(parRanges))
colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## add unique hash identifier
## (at the moment don't use "a0" type ensembleID, because MetaWards
## parses to dates)
design$output <- ensembleIDGen(ensembleID = "Ens0", nrow(design))
design$repeats <- 2

## convert to input space
input <- convertDesignToInput(design, parRanges, "zero_one")

## convert input to disease
disease <- convertInputToDisease(input, C)

## write to external files
dir.create("inputs", showWarnings = FALSE)

## write text file for MetaWards
write.table(disease, "inputs/disease.dat", row.names = FALSE, sep = " ", quote = FALSE)

## save inputs for data for post-simulation runs
saveRDS(design, "inputs/design.rds")
saveRDS(parRanges, "inputs/parRanges.rds")
