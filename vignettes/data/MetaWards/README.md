# MetaWards model construct

This directory contains the information necessary to hopefully run the
specific version of the MetaWards model specified. The details of the model
and explanations of the code can be found in the accompanying vignette.

This assumes you have **MetaWards version 1.2.0** installed; if you don't then
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

Once the design has been run, the `extractOutputsSQL.R` file provides some code
to extract some summary measures from the outputs and build an SQL database holding
the summary measures for each design point / replicate / ward combination. This is
stored in `raw_outputs/uberStages.db.bz2`.

This SQL database can be queried through any SQLite
client. The `dplyr` package (or more specifically the `dbplyr` package) provides
some useful R tools for querying SQL databases using `tidyverse`-type notation. Some
examples are in the `extractOutputs.R` file, and more details can be found
[here](https://cran.r-project.org/web/packages/dbplyr/vignettes/dbplyr.html).

There is subsequently a file called `extractOutputs.R` that includes some R code 
giving examples about how to manipulate and extract information from `uberStages.db`.
Note that in order to do this you will have to uncompress the `uberStages.db.bz2` file.
On Linux this can be done on the command line, e.g.

```
bzip2 -dkf uberStages.db.bz2
```

For large designs the uncompressed file will be **huge**, and it will take a while to 
decompress. Whether you choose to store the decompressed file permanently will depend
on how much storage you have!

## Animations

In the `images` folder there are some R scripts to produce animations. The `plotAnimation.R` script accesses weekly counts from the `uberStages.db` database. It takes three inputs:

* ID (Ensemble ID, e.g. "Ens0000")
* REP (Replicate number)
* VAR (Variable you wish to animate---must be in `uberStages.db`)

So the command:

```
R CMD BATCH --no-restore --no-save --slave "--args Ens0000 1 Hprev" plotAnimation.R
```

will produce an animation of the weekly `Hprev` values for replicate 1 of design point `Ens0000`. **Note: it helps to have set an index on the `output` column of `uberStages.db`**---see the comments in `plotAnimation.R` for more details.

Alternatively, the `plotAnimation_stages.R` script accesses daily counts from an individual the `stages.db.bz2` file. It takes three inputs:

* ID (Ensemble ID, e.g. "Ens0000")
* REP (Replicate number)
* VAR (Variable you wish to animate---must be in `stages.db`)

So the command:

```
R CMD BATCH --no-restore --no-save --slave "--args Ens0000 1 H" plotAnimation_stages.R
```

will produce an animation of the daily `H` values for replicate 1 of design point `Ens0000`. This is quicker due to using base R plotting, rather than `gganimate`.

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
