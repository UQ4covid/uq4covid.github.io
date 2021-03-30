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
and the various design files edited at will as required.

Similarly, generating the vignette also creates a `.zip` file called `metawardsCatalyst.zip`
that contains the same, but with additional code for running on Catalyst and post-processing
on JASMIN. **This will not work unless you have all the necessary datasets described above,
some of which are not publicly available.** Once the design has been included, this can be sent 
to Christopher, who will run it using the `catalystJobscript.sh` script and copy the results 
to JASMIN in the `covid19` workspace. 

## JASMIN post-processing

### Logging on to JASMIN

Assuming you're set up on JASMIN, and have access to the GWS for the `covid19`
project, then if you're not on the University network, then something like the following 
should work (of course, replacing relevant login details and paths):

```
eval $(ssh-agent -s)
ssh-add PATH_TO_PRIVATE_KEY
ssh -A USERNAME@login2.jasmin.ac.uk
```

For more details about how to set-up an account on JASMIN, see 
[https://help.jasmin.ac.uk/article/189-get-started-with-jasmin](https://help.jasmin.ac.uk/article/189-get-started-with-jasmin).

### Install necessary R libraries

**Note**: in order to run R scripts on JASMIN you will have to ensure that you install the necessary
R libraries in your personal libraries. Once you're on the login node as described above, `ssh` 
into one of the scientific processing nodes e.g.

```
ssh USERNAME@sci1.jasmin.ac.uk
```

Then load `jaspy`, which gives access to R:

```
module load jaspy
```

Now you should load R and then install whatever packages you need in the usual way. You will need at
least `tidyverse` and `RSQLite` for the code below to work (e.g. `install.packages("tidyverse", "RSQLite")`).

### Copying files to JASMIN

Note that all runs processed on Catalyst will be copied into `/gws/nopw/j04/covid19/catalyst`
after they have been run, in which case you can skip this section.

Assuming you're set up on JASMIN, and have access to the GWS for the `covid19`
project. Then you need to login to a transfer server to copy the files across.

Once you're connected to the login node, you can `ssh` into the transfer node as:

```
ssh USERNAME@xfer1.jasmin.ac.uk
```

From here, change directory to the `covid19` GWS workspace:

```
cd /gws/nopw/j04/covid19/
```

**NOTE**: you can only copy files **to** JASMIN if you are connected to a network that
allows access to JASMIN (e.g. the University of Exeter network). This won't work
over a VPN. You can connect from JASMIN to another external IP via `scp` or `rsync`, e.g.

```
scp USER@REMOTEIP:PATHTOFOLDER .
```

Once this is done, disconnect from the `xfer*` node.

```
exit
```

### Extract all `.zip` files

**Note**: the extraction code below only has to happen once, and for any runs from Catalyst
I will aim to do this soon after the runs are transferred, hence you can skip this section.

Firstly, from the login node, log in to one of the scientific processing nodes e.g.

```
ssh USERNAME@sci1.jasmin.ac.uk
```

Then change directory to the relevant folder containing the JASMIN code e.g.

```
cd /gws/nopw/j04/covid19/FOLDER/JASMINsetup
```

where `FOLDER` is replaced with the correct folder path e.g. `catalyst/wave0`.

Now load the `jaspy` module:

```
module load jaspy
```

Now edit the `setupSLURM.R` script and change `filedir` to create a directory in the `public`
folder where the outputs will be stored. This will be made a
sub-directory of `/gws/nopw/j04/covid19/public`. You will also need
to change the `startdate` to match the start of the MetaWards simulation, and 
`ndays` to reflect the number of days that the simulation is run over e.g.

```
## set directory to save outputs to and startdate
filedir <- "wave0"
startdate <- "09/02/2020"
ndays <- 41
```

Then run the script:

```
R CMD BATCH --no-restore --no-save --slave setupSLURM.R
```

This will create a file called `submit_job.sbatch` that we can submit to SLURM via:

```
sbatch submit_job.sbatch
```

This will run once the scheduler allows. If you want to change any of the settings (like
wall time etc.), then either edit the `submit_job.sbatch` file directly,
or alter the `submit_job_template.sbatch` template file and then re-run `setupSLURM.R` as
above.

After all of these jobs have completed, the unzipped SQLite databases can be found on the
public repository.

### Producing summary tables on JASMIN

**Note**: We have provided a template folder `JASMINsummary` that contains code that
can be amended to produce different summary measures as required. In the discussion
below we will extract weekly average hospital prevalences and weekly
cumulative deaths for all wards and weeks.

If you wish to write your own post-processing code, it is a good idea to copy
the `JASMINsummary` folder and rename it, and then make changes to the code
as detailed below e.g.

```
cp -r JASMINsummary JASMINsummary_new
```

**Note:** for the code to work you must not move your copy of `JASMINsummary` to any other
location, since various relative links are used throughout for ease.

The discussion below will assume that you are just running the default 
post-processing code.

Firstly, from the login node, log in to e.g `sci1.jasmin.ac.uk` as before:

```
ssh USERNAME@sci1.jasmin.ac.uk
```

Then change directory to the relevant folder containing the JASMIN code e.g.

```
cd /gws/nopw/j04/covid19/FOLDER/JASMINsummary_new
```

where `FOLDER` is replaced with the correct folder path e.g. `catalyst/wave0`.

Now load the `jaspy` module:

```
module load jaspy
```

Now edit the `setupSLURM.R` script and set a unique identifier that will be appended 
to your summary runs (see below).

**Note**: the unique user ID is important to stop your outputs from overwriting 
anyone else's.

Hence amend the line below as appropriate:

```
## set unique ID
id <- "user"
```

If you wanted to run the script for a subset of inputs, you could filter
the `inputs` data frame accordingly in this file directly after reading it in.

Once you're happy, you can run the script:

```
R CMD BATCH --no-restore --no-save --slave setupSLURM.R
```

This will create a file called `submit_job.sbatch` that we can use to submit jobs
to SLURM.

Before we do that, the script file `createSum.R` contains the code necessary
to extract the necessary summary measures. You should be able to edit
lines 85--97 accordingly to generate different summary measures as required.

Once you're happy, the jobs can be submitted via:

```
sbatch submit_job.sbatch
```

This will run once the scheduler allows. If you want to change any of the settings (like
wall time etc.), then either edit the `submit_job.sbatch` file directly,
or alter the `submit_job_template.sbatch` template file and then re-run `setupSLURM.R` as
above.

This creates a data frame, saved as a `output_ID.rds` file in each folder on the public repo
(where `ID` is replaced with the unique identifier from your `setupSLURM.R` script). Once 
all these jobs have completed, you may want to run a script to collate these results together
in one data frame or file for ease of downloading / querying. To this end there are two
auxiliary files: `collateSum.R` and `collateSumSQL.R`. 

**Important**: the former collates all the summaries into a single `.csv` file that is stored 
in the main repo e.g. 
`/gws/nopw/j04/covid19/public/wave0/raw_outputs/summaries_ID.csv` (where `ID` is
the explicit unique identifier used in the `setupSLURM.R` script). Please note, only
do this if you have created summaries at some higher spatial resolution (e.g. trusts / LADs
and not wards), otherwise the file will be too large to create on JASMIN and make
querying more difficult. 

For ward-level summaries, it is better to use the `collateSumSQL.R` script, which
collates all the summaries into a single compressed SQLite database 
that is stored in the main repo e.g. 
`/gws/nopw/j04/covid19/public/wave0/raw_outputs/summaries_ID.db` (where `ID` is
the explicit unique identifier used in the `setupSLURM.R` script) e.g.

```
R CMD BATCH --no-restore --no-save --slave collateSumSQL.R
```

Once this has been checked, you can run the `cleanup.R` script to remove all the
intermediate `output_ID.rds` data frames and clean up the repo. It's a good idea to wait 
until you're sure the collate function has run properly before cleaning up these files.
You can easily amend the cleanup file to leave the intermediate data frames as required.

```
R CMD BATCH --no-restore --no-save --slave cleanup.R
```

## Querying files from external sources

All relevant files should be accessible on a server than can be directly accessed
at [https://gws-access.jasmin.ac.uk/public/covid19/](https://gws-access.jasmin.ac.uk/public/covid19/). You cannot query an SQLite database from a server like this, you can only
download files. 

