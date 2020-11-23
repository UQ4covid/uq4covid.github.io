## load libraries
library(lhs)
library(dplyr)
library(tidyr)

## source dataTools
source("R_tools/dataTools.R")

## set up parameter ranges
parRanges <- data.frame(
    parameter = c("r_zero", "incubation_time", "infectious_time", "hospital_time",
                  "critical_time", "lock_1_restrict", "lock_2_release",
                  "pEA", "pIH", "pIRprime", "pHC", "pHRprime", "pCR", 
                  "GP_A", "GP_H", "GP_C"),
    lower = c(2, 4, 2, 4, 4, rep(0, 11)),
    upper = c(4, 6, 4, 12, 12, rep(1, 11)),
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
