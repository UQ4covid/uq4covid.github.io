# MetaWards model construct

This directory contains the information necessary to hopefully run the
specific version of the MetaWards model specified. The details of the model
and explanations of the code can be found in the accompanying vignette.

This assumes you have MetaWards version 1.1.0 installed; if you don't then
you can find instructions [here](https://metawards.org/install.html). It also
assumes that you have cloned the MetaWardsData as detailed 
[here](https://metawards.org/model_data.html).

The current code is written to run on Linux, but it should also run happily 
on a Mac (not tested yet). It should also run on Windows, but you will have to convert
`runscript.sh` to a Windows batch file (Pull Requests gratefully received).

In `runscript.sh` you will need to amend the line:

```
export METAWARDSDATA=$HOME/Documents/covid/MetaWardsData
```

to point to your installation on the MetaWardsData repository.

The R script `convertDesign.R` contains a simple example of sampling
from a $(0, 1)$ LHS sampler, and then converting the design into
the correct format for MetaWards, which is written into the `inputs/disease.dat`
file. This uses `convertDesigntoInput()` and `convertInputToDisease()` tools
that can be found in the `R_tools/dataTools.R` script, and are described in the
vignette. (These tools also accept designs in $(-1, 1)$ spaces.)

Once the design has been run, the `extractOutputs.R` file provides some code
to extract some summary measures from the outputs. At the current time this relies
on some OS-specific code, in order to unzip the compressed SQL database files. Hopefully
this will be sorted with a better solution in the future, but it gives you some
ideas about what can be done. This is also documented in the vignette.

The outputs are stored in an SQL database, which can be queried through any SQLite
client. The `dplyr` package (or more specifically the `dbplyr` package) provides
some useful R tools for querying SQL databases using `tidyverse`-type notation. Some
examples are in the `extractOutputs.R` file, and more details can be found
[here](https://cran.r-project.org/web/packages/dbplyr/vignettes/dbplyr.html).

## Possible extensions / to-do

* Lockdown cut-off for distance travelled.
* Amend lockdown iterator to model weekdays and weekends during lockdown.
* Superspreaders / supershedders?
* Possible additional hospital workers class?
* Change names of outputs to something easier to understand.
* Sort out how to generically unzip files rather than using `system()` (hopefully Chris' extractor will solve this).
* Need some checks of inputs in R tools.
* Perhaps come up with a better way to store the data (maybe only store days where some events have changed,
  and then post-process to fill in the gaps where necessary).
