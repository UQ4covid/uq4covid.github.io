The script files and data sources used here are to generate plausible ranges for the age-specific transition probabilities.

The data (`hospitalisationImperial.csv`) is a copy of Table 3 from [Verity et al., *Lancet Infectious Diseases*, 2020](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext). The data (`IFRImperial.csv`) is a copy of Table 1 from [Verity et al., *Lancet Infectious Diseases*, 2020](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext). 

The estimates of the infection fatality risk and the probability of hospitalisation given infection from the CDC (`CasesHospDeathInHospForUS_Update.xlsx`) are  obtained in the following way.

1. The cumulative number of cases and deaths for each age group by 26/02/2021 were downloaded from [https://covid.cdc.gov/covid-data-tracker/#demographics](https://covid.cdc.gov/covid-data-tracker/#demographics).
2. We use the cumulative hospitalisation proportion (per 100,000 population) on 20/02/2021 (available from [https://gis.cdc.gov/grasp/COVIDNet/COVID19_3.html](https://gis.cdc.gov/grasp/COVIDNet/COVID19_3.html)), and the US population by age (available from [https://www.statista.com/statistics/241488/population-of-the-us-by-sex-and-age/](https://www.statista.com/statistics/241488/population-of-the-us-by-sex-and-age/)) to estimate the number of hospitalisations for the entire US population.
3. We obtain the number of in hospital deaths by age group from [https://data.cdc.gov/dataset/NVSS-Provisional-COVID-19-Deaths-by-Place-of-Death/4va6-ph5s/data](https://data.cdc.gov/dataset/NVSS-Provisional-COVID-19-Deaths-by-Place-of-Death/4va6-ph5s/data), by filtering the data for the entire US using the category "Place of death" / "Healthcare setting, inpatient". 

To run the analysis, set the working directory to this folder, and then run the `cleanData.R` script file to clean and combine data for the subsequent analyses.

Running the `sims.R` script file reruns our plausible range code, and produces a file called `pathways.rds` that contains a finite mixture model that can be used to generate design points for these parameters. By default this is copied to the `../../inputs/` folder.

Then running the `deriveLimits.R` file produces thresholds for the mixture densities that can be used to reject
points in later waves. The object `pathThresh.rds` is added to the `../../inputs/` folder.