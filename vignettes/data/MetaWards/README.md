# MetaWards model construct

This directory contains the information necessary to hopefully run the
specific version of the MetaWards model specified. The details of the model
and explanations of the code can be found in the accompanying vignette.

This assumes you have **MetaWards version 1.5.1** installed; if you don't then
you can find instructions [here](https://metawards.org/install.html). It also
assumes that you have cloned the MetaWardsData as detailed 
[here](https://metawards.org/model_data.html).

## Setting up data

Some data generation must be done in advance as documented below.

The folder `data/ward2011ToWard2019Mapping` folder contains 
the relevant shapefiles, lookup tables and code to convert
the 2011 ward-level commuter data to the 2019 wards and LADs. There is also
code to generate the relevant data for use in the MetaWards model runs.

The `data/hospitalStays` folder contains information about
how to inform the length of hospital stays from CHESS data. (See `data/hospitalStays/README.md`.)

The input data `data/populationByAge/Pop by CTRY.csv` is from ONS census and contains 
population counts in different age-classes in different countries. (See `data/populationByAge/README.md`.)

The file `README.md` in the `data/contactData` folder contains information about
how to download and generate the contact matrices used in the model.

The `data/hospitalCatchments` folder contains scripts provided by Rob Challen for mapping
NHS Trusts to 2019 wards.

The `data/intervals` folder contains relevant screenshots from [Challen et al. (2020)](https://www.medrxiv.org/content/10.1101/2020.11.17.20231548v2).

The `data/pathways` folder contains data from [Verity et al. (2020)](https://www.thelancet.com/journals/laninf/article/PIIS1473-3099(20)30243-7/fulltext)
and the CDC website, along with code to generate plausible ranges for various
transition probability parameters.

The `data/seedDeaths` folder contains code to generate seeding times and probabilities.
This relies on death data we are not able to share publicly.

## Setting up design

The R script `convertDesign.R` contains a ***simple example*** of sampling
from a $(0, 1)$ LHS sampler, and then converting the design into
the correct format for MetaWards, which is written into the `inputs/disease.dat`
file. This uses `convertDesigntoInput()` and `convertInputToDisease()` tools
that can be found in the `R_tools/dataTools.R` script, and are described in the
vignette. (These tools also accept designs in $(-1, 1)$ spaces as long as you
set the `scale` argument to `convertDesignToInput()`.)

Other designs can be added using this code as a guide.

### Tests

There are some tests in the `tests/` folder. The `bash` script file `testScript.sh`
gives details. These run tests in a single population, and compare outputs to a deterministic
model (to check that the finalsizes match those that are obtained from the NGM matrix
using the approach of [Andreasen (2011)](https://link.springer.com/article/10.1007%2Fs11538-010-9623-3)). 
They also compare the MetaWards runs to those obtained form an equivalent single
population model coded in R in order to check various aspects of the implementation. 
Paths to the `MetaWardsData` folder will have to be updated on different machines.

## Running the model on Catalyst

When the vignette is generated, this produces a `.zip` file called `metawards.zip` that
contains all the code necessary to run the model on an individual machine. This can be unzipped
and the varius design files edited at will as required.

Similarly, generating the vignette also creates a .zip. file called `metawardsCatalyst.zip`
that contains the same, but with additional code for running on Catalyst and post-processing
on JASMIN. Once the design has been included, this can be sent to Christopher, who will run it 
using the `catalystJobscript.sh` script and copy the results to JASMIN in the `covid19` workspace. 

## JASMIN post-processing

### Copying files to JASMIN

Assuming you're set up on JASMIN, and have access to the GWS for the `covid19`
project. Then you need to login to a transfer server to copy the files across.
If you're not on the University network, then something like the following should 
work (of course, replacing relevant login details and paths):

```
eval $(ssh-agent -s)
ssh-add PATH_TO_PRIVATE_KEY
ssh -A USERNAME@login2.jasmin.ac.uk
```

Once you're connected to the login node, you can `ssh` into the transfer node as:

```
ssh USERNAME@xfer1.jasmin.ac.uk
```

From here, change directory to the `covid19` GWS workspace:

```
cd /gws/nopw/j04/covid19/
```

Transfer of runs to this folder can be done via `scp` or `rsync`, e.g.

```
scp USER@REMOTEIP:PATHTOFOLDER .
```

Once this is done, disconnect from the `xfer*` node.

```
exit
```

### Create script

Firstly, from the login node, log in to `sci1.jasmin.ac.uk` e.g.

```
ssh USERNAME@sci1.jasmin.ac.uk
```

Then change directory to the relevant folder containing the JASMIN code e.g.

```
cd /gws/nopw/j04/covid19/FOLDER/JASMIN
```

Now load the `jaspy` module:

```
module load jaspy
```

Now edit the `setupSLURM.R` script and change `filedir` to point to the public
directory that you want the files saved into (with trailing `/`) e.g.

```
filedir <- "/gws/nopw/j04/covid19/public/wave0/"
```

This must be a sub-directory of `/gws/nopw/j04/covid19/public`. Then run the script:

```
R CMD BATCH --no-restore --no-save --slave setupSLURM.R
```

This will create a file called `submit_job.sbatch` that we can submit to SLURM. 

### Extracting outputs and producing summary tables on JASMIN

This next step can be done in parallel, and can be run by submitting a batch
job script via SLURM:

```
sbatch submit_job.sbatch
```

This will run once the scheduler allows. If you want to change any of the settings (like
wall time etc.), then either edit the `submit_job.sbatch` file directly,
or alter the `submit_job_template.sbatch` template file and then re-run `setupSLURM.R` as
above.

## Querying files from external sources

All relevant files should be accessible on a server than can be directly accessed
at [https://gws-access.jasmin.ac.uk/public/covid19/](https://gws-access.jasmin.ac.uk/public/covid19/). You cannot query an SQLite database from a server like this, you can only
download files. Thus the `age*.db` databases in each sub-directory contain the
raw outputs, but the summary measures are in the `weeksums_*.csv` files. These hold 
weekly average hospital prevalence (`Hprev`), along with the number of hospital deaths
(`Hdeaths`) and other deaths (`Cdeaths`) for each week / ward combination.

