The script `seedDeaths.R` produces an input table for seeding the model. 
It needs to be run AFTER `../hospitalCatchments/hospitalCatchments.R` has been run.

The data `ltla_2021-05-18.csv` are from:

[https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=cumDeathsByDeathDate&format=csv](https://api.coronavirus.data.gov.uk/v2/data?areaType=ltla&metric=cumDeathsByDeathDate&format=csv)

LAD to region lookup table `Local_Authority_District_to_Region_(April_2019)_Lookup_in_England.csv` can be found here:

[https://geoportal.statistics.gov.uk/datasets/local-authority-district-to-region-april-2019-lookup-in-england](https://geoportal.statistics.gov.uk/datasets/local-authority-district-to-region-april-2019-lookup-in-england)