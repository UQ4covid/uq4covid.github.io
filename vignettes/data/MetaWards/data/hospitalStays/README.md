The script files and data sources used here are to inform hospital length distributions by age.

Data (`admissionToOutcomeLineList.csv`) from Rob Challen.

Set working directory to this folder, and then the `hospitalStays.R` script file produces a file called `hospStays.rds` that contains a finite mixture model that can be used to generate design points for these parameters. By default this is copied to the `../../inputs/` folder.
