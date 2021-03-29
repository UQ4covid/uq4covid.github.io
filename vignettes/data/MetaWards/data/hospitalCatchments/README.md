The `WD19toNHSTrust` and `WD19toNHSTrustSupply` shapefiles are generated from Rob Challen using a novel method for generating NHS catchment areas. This also outputs the `WD19toNHSTrustCode.csv` lookup file.

Details on Rob's method can be found here:

[https://terminological.github.io/arear/articles/catchment-areas.pdf](https://terminological.github.io/arear/articles/catchment-areas.pdf)

Full code and documentation for the method can be found:

[https://github.com/terminological/arear](https://github.com/terminological/arear)

[https://terminological.github.io/arear/](https://terminological.github.io/arear/)

The shapefiles provided must be unzipped before the script described below is run. The `../Ward2011toWard19Mapping/create2011to2019.R` script file must also be run, and the results added to the `MetaWardsData` repository on the local machine before the steps below.

The path to the MetaWardsData directory might need to be amended in the `hospitalCatchments.R` script file, but running this file will produce a `../../inputs/trust19Lookup.csv` file that provides a way to map the 2019 wards to NHS Trusts for calibration / mapping purposes.

