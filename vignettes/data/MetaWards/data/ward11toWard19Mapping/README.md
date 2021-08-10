The script files and data sources used here allow for the 2011 ward-level commuting data to
be proportionately redistributed to the 2019 wards.

2011 Ward Lookup from the `MetaWardsData` repository. Need to add correct path to local repository
in the `create2011to2019.R` source file before running.

Ward19 Lookup (`Ward_to_Local_Authority_District_(December_2019)_Lookup_in_the_United_Kingdom.csv`) from:

(https://geoportal.statistics.gov.uk/datasets/ward-to-local-authority-district-december-2019-lookup-in-the-united-kingdom)[https://geoportal.statistics.gov.uk/datasets/ward-to-local-authority-district-december-2019-lookup-in-the-united-kingdom]

Shapefiles for 2011 Wards and 2019 Wards from:

(https://geoportal.statistics.gov.uk/datasets/wards-december-2011-boundaries-ew-bfc)[https://geoportal.statistics.gov.uk/datasets/wards-december-2011-boundaries-ew-bfc]
(https://geoportal.statistics.gov.uk/datasets/wards-december-2019-boundaries-ew-bfc)[https://geoportal.statistics.gov.uk/datasets/wards-december-2019-boundaries-ew-bfc]

Shapefiles for 2011 LADs and 2019 LADs from:

(https://geoportal.statistics.gov.uk/datasets/local-authority-districts-december-2011-boundaries-ew-bfc)[https://geoportal.statistics.gov.uk/datasets/local-authority-districts-december-2011-boundaries-ew-bfc]
(https://geoportal.statistics.gov.uk/datasets/local-authority-districts-december-2019-boundaries-uk-bfc)[https://geoportal.statistics.gov.uk/datasets/local-authority-districts-december-2019-boundaries-uk-bfc]

The ZIP files must be unzipped and all data files (except the 2011 Ward Lookup) must be
in the `ward11toWard19Mapping` folder, and the working directory needs to also be set to this folder
before running `create2011to2019.R`. This creates a new folder called `2011to2019Data` that 
**must** be added to the `MetaWardsData/model_data/` repo. In addition, the ward lookup table 
is automatically copied to `../../inputs/Ward19_Lookup.csv`.

The script `checks.R` runs a few manual comparisons as a sanity check.
