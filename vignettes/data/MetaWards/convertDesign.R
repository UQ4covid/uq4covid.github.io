## load libraries
library(lhs)
library(dplyr)
library(tidyr)

## source dataTools
source("dataTools.R")

## set up parameter ranges
parRanges <- data.frame(
    parameter = c("r_zero", "incubation_time", "infectious_time", "hospital_time",
                  "critical_time", "lock_1_restrict", "lock_2_release",
                  "pEA", "pIH", "pIRprime", "pHC", "pHRprime", "pCR", 
                  "GP_A", "GP_H", "GP_C"),
    lower = c(2.5, 4, 2, 4, 4, rep(0, 11)),
    upper = c(4, 6, 4, 12, 12, rep(1, 11)),
    stringsAsFactors = FALSE
)
print("NEED TO GET RANGES FIGURED OUT, THIS IS A TEST")

## generate LHS design
design <- randomLHS(5, nrow(parRanges))
colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## convert to input space
input <- convertDesignToInput(design, parRanges, "zero_one")

## convert input to disease
disease <- convertInputToDisease(input, 2)

## write to external files
dir.create("inputs")
write.table(disease, "inputs/disease.dat", row.names = FALSE, sep = " ", quote = FALSE)

