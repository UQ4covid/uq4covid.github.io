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
  lower = c(2.5, 4, 2, 4, 4, rep(0, 11)),
  upper = c(4, 6, 4, 12, 12, rep(1, 11)),
  stringsAsFactors = FALSE
)

## read in Danny's design on [0,1]
source("kExtendedLHCs.R")
Rep_Ens_Size <- 40 #? 15 parameters so if had 180 runs, can leave out 30 at a time and still have 10p
HowManyCubes <- 5
New_cube <- MakeRankExtensionCubes(n=Rep_Ens_Size, m=16, 
                                   k=HowManyCubes, w=0.2, FAC_t=0.5)
DannyDesign200 <- NewExtendingLHS(New_cube)
saveRDS(DannyDesign200, "DannyDesign200.rds")

design <- readRDS("DannyDesign200.rds")

colnames(design) <- parRanges$parameter
design <- as_tibble(design)

## add unique hash identifier
## (at the moment don't use "a0" type ensembleID, because MetaWards
## parses to dates)
design$output <- ensembleIDGen(ensembleID = "Ens1", nrow(design))
design$repeats <- 20

## convert to input space
input <- convertDesignToInput(design, parRanges, "zero_one")

## convert input to disease
disease <- convertInputToDisease(input)

## write to external files
dir.create("inputs", showWarnings = FALSE)

## write text file for MetaWards
write.table(disease, "inputs/disease.dat", row.names = FALSE, sep = " ", quote = FALSE)

## save inputs for data for post-simulation runs
saveRDS(design, "inputs/design.rds")
saveRDS(parRanges, "inputs/parRanges.rds")
