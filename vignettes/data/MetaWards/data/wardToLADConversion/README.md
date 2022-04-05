Shapefiles for 2019 LADs from:

[https://geoportal.statistics.gov.uk/datasets/local-authority-districts-december-2019-boundaries-uk-buc/explore](https://geoportal.statistics.gov.uk/datasets/local-authority-districts-december-2019-boundaries-uk-buc/explore)

This takes the data from the `MetaWardsData/model_data/2011to2019Data` repo and converts to LADs.

The shapefile above must be unzipped and located in the `wardtoLADConversion` folder, and the working
directory needs to also be set to this folder before running `wardToLADConversion.R`. 
This creates a new folder called `2019LADData` that **must** be added to the 
`MetaWardsData/model_data/` repo. In addition, the ward lookup table is automatically copied 
to `../../inputs/Ward19_Lookup.csv`.
