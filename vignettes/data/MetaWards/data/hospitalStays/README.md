The script files and data sources used here are to inform hospital length distributions by age.

Data are from the CHESS study which we are not able to share publicly, so we make the
code but not the data available.

'SGSS and CHESS data', NHS Digital. [https://digital.nhs.uk/about-nhs-digital/corporate-information-and-documents/directions-and-data-provision-notices/data-provision-notices-dpns/sgss-and-chess-data](https://digital.nhs.uk/about-nhs-digital/corporate-information-and-documents/directions-and-data-provision-notices/data-provision-notices-dpns/sgss-and-chess-data) (accessed Dec. 18, 2020).

To run the analysis, set working directory to this folder, and then the `hospitalStays.R` script file produces a file called `hospStays.rds` that contains a finite mixture model that can be used to generate design points for these parameters. By default this is copied to the `../../inputs/` folder.
