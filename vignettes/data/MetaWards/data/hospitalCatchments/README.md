The `WD19toNHSTrust` and `WD19toNHSTrustSupply` shapefiles are generated from Rob Challen using his NHS catchment area method. This also outputs the `WD19toNHSTrustCode.csv` lookup file.

The shapefiles must be unzipped before the script described below is run. The `../Ward2011toWard19Mapping/create2011to2019.R` script file must also be run, and the results added to the `MetaWardsData` repository on the local machine before the steps below.

The path to the MetaWardsData directory might need to be amended in the `hospitalCatchments.R` script file, but running this file will produce a `trust19Lookup.csv` file that provides a way to map the 2019 wards to NHS Trusts for calibration / mapping purposes.

