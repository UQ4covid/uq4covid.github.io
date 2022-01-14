This folder contains a simulation study to test the PF
approach to calibrate models and perform DA simultaneously.

The files `disease.dat` and `inputs.rds` contain the design
we have been using for the MW model. Since we are running
without age-structure here, the transmission parameter
needs adjusting for the lack of age-structure.

Firstly run `stochModel.R` to simulate data using one
of the MW design points. Outputs are placed in the `outputs`
folder (which is created if it is not present).

Then the file `analysis.R` gives some commands for reading
in the data and the current design (with some conversion to
the `nu` parameter to account for the lack-of-age-structure).

Some code is shown for running the model and PF and generating
log-likelihood estimates for the design, as well as producing
some plots of the particle densities (unweighted, so only an 
approximation here) against the data for comparison, if
you wish to do this. The PF code is in the `PF.R` file, the 
truncated Skellam sampler is in the `trSkellam.R` file
and the next-generation matrix function to convert `nu`
is in `NGM.R` (you're probably less interested in this).
