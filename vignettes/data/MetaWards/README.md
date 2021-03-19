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

The file `README.md` in the `data/ward2011ToWard2019Mapping` folder contains information about
where to download the relevant shapefiles and lookup tables to convert
the 2011 ward-level commuter data to the 2019 wards and LADs. There is also
code to generate the relevant data for use in the MetaWards model runs.

The file `README.md` in the `data/hospitalStays` folder contains information about
how to inform the length of hospital stays from line list data provided by Rob Challen.

The input data `data/populationByAge/Pop by CTRY.csv` is from Rob Challen and contains 
population counts in different age-classes in different countries.

The file `README.md` in the `data/contactData` folder contains information about
how to download and generate the contact matrix data used in the model.

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
that contains the same, but with additional code for runnign on Catalyst and post-processing
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

### Create scripts

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

Now edit the `setupLOTUS.R` script and change `filedir` to point to the public
directory that you want the files saved into (with trailing `/`) e.g.

```
filedir <- "/gws/nopw/j04/covid19/public/wave0/"
```

This must be a sub-directory of `/gws/nopw/j04/covid19/public`. Then run the script:

```
R CMD BATCH --no-restore --no-save --slave setupSLURM.R
```

This will create two files called `submit_job.sbatch` and `submit_quantile.sbatch`
that we can submit to SLURM. Edit the `createSum.sh` and `createQuantiles.sh` files 
and change the lines:

```
filedir="/gws/nopw/j04/covid19/public/wave0/"
```

to point to the correct output folder that we wish other people to access (with
trailing `/` as before).

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

### Extracting quantiles

Once the above has been done, submit the following job file via SLURM:

```
sbatch submit_quantile.sbatch
```

This will run once the scheduler allows. If you want to change any of the settings (like
wall time etc.), then either edit the `submit_quantile.sbatch` file directly,
or alter the `submit_quantile_template.sbatch` template file and then re-run `setupSLURM.R` as
above.


## Querying files from external sources

All relevant files should be accessible on a server than can be directly accessed
at [https://gws-access.jasmin.ac.uk/public/covid19/](https://gws-access.jasmin.ac.uk/public/covid19/). You cannot query an SQLite database from a server like this, you can only
download files. Thus the `age*.db` databases in each sub-directory contain the
raw outputs, but the summary measures are in the `weeksums.csv` files. These hold 
weekly average `Hprev` and total `Deaths` for each week / ward combination 
for every week since just before the first lockdown. 

In the baseline directory for each design point (e.g. `Ens0000` and not `Ens0000x001` etc.)
there are two files: `Hprev.csv` and `Deaths.csv` that contain the
quantiles of the different outputs over the replicates for each design point.

The R script `downloadQuantiles.R` provides some parallel code that can be run from 
any user machine to download these and concatenate the outputs across the design
points accordingly. Just change the lines:

```
filedir <- "https://gws-access.jasmin.ac.uk/public/covid19/wave0/"
week <- 12
output <- "Hprev"
ncores <- 24
```

accordingly. The first line needs to point to the correct directory on the server 
(with trailing `/`). The `week` object can either be a vector of weeks to extract,
or `NA`, in which case it will extract all weeks. The `output` object is the name
of the output to extract (it must be one of: `Hprev` or `Deaths`). Finally,
the `ncores` object gives the number of cores to use to read entries in parallel
(Windows user will have to set this to be 1).

If you want to extract for a subset of design points, then swap `design$output` in 
the `mclapply()` function for a vector of design hashes.

