# MetaWards model construct

This directory contains the information necessary to hopefully run the
specific version of the MetaWards model specified. The details of the model
and explanations of the code can be found in the accompanying vignette.

This assumes you have **MetaWards version 1.2.0** installed; if you don't then
you can find instructions [here](https://metawards.org/install.html). It also
assumes that you have cloned the MetaWardsData as detailed 
[here](https://metawards.org/model_data.html).

## Setting up data

The file `README.md` in the `data/tierData` folder contains information about
where to download the relevant shapefiles and lookup tables to convert
the 2011 ward-level commuter data to the 2019 wards and LADs. There is also
code to generate the relevant data for use in the MetaWards model runs.

The file `README.md` in the `data/hospitalStays` folder contains information about
how to inform the length of hospital stays from line list data.

## Setting up design

The R script `convertDesign.R` contains a simple example of sampling
from a $(0, 1)$ LHS sampler, and then converting the design into
the correct format for MetaWards, which is written into the `inputs/disease.dat`
file. This uses `convertDesigntoInput()` and `convertInputToDisease()` tools
that can be found in the `R_tools/dataTools.R` script, and are described in the
vignette. (These tools also accept designs in $(-1, 1)$ spaces as long as you
set the `scale` argument to `convertDesignToInput()`.)

## Running the model on Catalyst

Once the design has been generated and the `inputs/disease.dat` file specified, the
whole folder can be zipped and sent to Christopher, who will run it using the
`jobscript.sh` script. 

You can amend this script if you want to change anything (for example, the 
`--nsteps 177` default number of days to run the model for might want amending).

Once the model has been run, Christopher will transfer to the AWS server.

## Copying files to JASMIN

Assuming you're set up on JASMIN, and have access to the GWS for the `covid19`
project. Then you need to login to a transfer server to copy the files across.
If you're not on the University network, then something like the following should 
work (of course, replacing relevant login details and paths):

```
eval $(ssh-agent -s)
ssh-add PATH_TO_PRIVATE_KEY
ssh -A USERNAME@jasmin-login2.ceda.ac.uk
```

Once you're connected to the login node, you can `ssh` into the transfer node as:

```
ssh USERNAME@jasmin-xfer1.ceda.ac.uk
```

From here, change directory to the `covid19` GWS workspace:

```
cd /gws/nopw/j04/covid19/
```

Now you can transfer the relevant folder from the AWS machine across (making sure 
it's not overwriting anything by changing folder paths if necessary). Note that
the `aws_exeter.pem` key is already in this folder for ease.

```
scp -i ~/aws_exeter.pem ubuntu@35.178.206.202:/home/ubuntu/FOLDERNAME .
```

Once this is done, disconnect from the `xfer*` node.

```
exit
```

## Create scripts

Firstly, from the login node, log in to `sci1.ceda.ac.uk` e.g.

```
ssh USERNAME@sci1.ceda.ac.uk
```

Then change directory to the relevant folder:

```
cd /gws/nopw/j04/covid19/FOLDER
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
R CMD BATCH --no-restore --no-save --slave setupLOTUS.R
```

This will create two files called `submit_job.sbatch` and `submit_comb.sbatch`
that we can submit to LOTUS. Edit the `createSum.sh` and 
`createQuantiles.sh` files and change the lines:

```
filedir="/gws/nopw/j04/covid19/public/wave0/"
```

to point to the correct output folder that we wish other people to access (with
trailing `/` as before).

## Extracting outputs and producing summary tables on JASMIN

This next step can be done in parallel, and can be run by submitting a batch
job script via SLURM:

```
sbatch submit_job.sbatch
```

This will run once the scheduler allows. If you want to change any of the settings (like
wall time etc.), then either edit the `submit_job.sbatch` file directly,
or alter the `submit_job_template.sbatch` template file and then re-run `setupLOTUS.R` as
above.

## Extracting quantiles

Once the above has been done, submit the following job file via SLURM:

```
sbatch submit_quantile.sbatch
```

This will run once the scheduler allows. If you want to change any of the settings (like
wall time etc.), then either edit the `submit_quantile.sbatch` file directly,
or alter the `submit_quantile_template.sbatch` template file and then re-run `setupLOTUS.R` as
above.


## Querying files from external sources

All relevant files should be accessible on a server than can be directly accessed
at [https://gws-access.jasmin.ac.uk/public/covid19/](https://gws-access.jasmin.ac.uk/public/covid19/). You cannot query an SQLite database from a server like this, you can only
download files. Thus the `stages.db` databases in each sub-directory contain the
raw outputs, but the summary measures are in the `weeksums.csv` files. These hold 
weekly average `Hprev`, `Cprev` and total `Deaths` for each week / ward combination 
for every week since just before the first lockdown. 

In the baseline directory for each design point (e.g. `Ens0000` and not `Ens0000x001` etc.)
there are three files: `Hprev.csv`, `Cprev.csv` and `Deaths.csv` that contain the
quantiles of the different outputs over the replicates for each design point.

The R script `downloadQuantiles.R` provides some parallel code that can be run from 
any user machine to download these and concatentae the outputs across the design
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
of the output to extract (it must be one of: `Hprev`, `Cprev` or `Deaths`). Finally,
the `ncores` object gives the number of cores to use to read entries in parallel
(Windows user will have to set this to be 1).

If you want to extract for a subset of design points, then swap `design$output` in 
the `mclapply()` function for a vector of design hashes.


<!--## Animations-->

<!--In the `images` folder there are some R scripts to produce animations. The `plotAnimation.R` script accesses weekly counts from the `uberStages.db` database. It takes three inputs:-->

<!--* ID (Ensemble ID, e.g. "Ens0000")-->
<!--* REP (Replicate number)-->
<!--* VAR (Variable you wish to animate---must be in `uberStages.db`)-->

<!--So the command:-->

<!--```-->
<!--R CMD BATCH --no-restore --no-save --slave "--args Ens0000 1 Hprev" plotAnimation.R-->
<!--```-->

<!--will produce an animation of the weekly `Hprev` values for replicate 1 of design point `Ens0000`. **Note: it helps to have set an index on the `output` column of `uberStages.db`**---see the comments in `plotAnimation.R` for more details.-->

<!--Alternatively, the `plotAnimation_stages.R` script accesses daily counts from an individual the `stages.db.bz2` file. It takes three inputs:-->

<!--* ID (Ensemble ID, e.g. "Ens0000")-->
<!--* REP (Replicate number)-->
<!--* VAR (Variable you wish to animate---must be in `stages.db`)-->

<!--So the command:-->

<!--```-->
<!--R CMD BATCH --no-restore --no-save --slave "--args Ens0000 1 H" plotAnimation_stages.R-->
<!--```-->

<!--will produce an animation of the daily `H` values for replicate 1 of design point `Ens0000`. This is quicker due to using base R plotting, rather than `gganimate`.-->

<!--## Possible extensions / to-do-->

<!--* Lockdown cut-off for distance travelled.-->
<!--* Amend lockdown iterator to model weekdays and weekends during lockdown.-->
<!--* Superspreaders / supershedders?-->
<!--* Possible additional hospital workers class?-->
<!--* Change names of outputs to something easier to understand.-->
<!--* Sort out how to generically unzip files rather than using `system()` (hopefully Chris' extractor will solve this).-->
<!--* Need some checks of inputs in R tools.-->
<!--* Perhaps come up with a better way to store the data (maybe only store days where some events have changed,-->
<!--  and then post-process to fill in the gaps where necessary).-->
<!--  -->
