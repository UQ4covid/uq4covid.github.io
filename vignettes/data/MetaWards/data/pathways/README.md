The script files and data sources used here are to generate plausible ranges for the age-specific transition probabilities.

The data (`hospitalisationImperial.csv`) is a copy of Table 3 from Verity et al., *Lancet Infectious Diseases*, 2020. The data (`IFRImperial.csv`) is a copy of Table 1 from Verity et al., *Lancet Infectious Diseases*, 2020. The data (`CasesHospDeathInHospForUS_Update.xlsx`) is derived from the CDC website and other sources (**TO COMPLETE**).

Set working directory to this folder, and then run the `cleanData.R` script file to clean and combine data for the subsequent analyses.

Running the `sims.R` script file reruns our plausible range code, and produces a file called `pathways.rds` that contains a finite mixture model that can be used to generate design points for these parameters.
