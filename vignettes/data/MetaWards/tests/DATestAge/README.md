This folder contains a simulation study to test the PF
approach to calibrate models and perform DA simultaneously.

The files `disease.dat` and `inputs.rds` contain the design
we have been using for the MW model. 

Firstly run `stochModel.R` to simulate data using one
of the MW design points. Outputs are placed in the `outputs`
folder (which is created if it is not present).

Then the file `analysis.R` gives some commands for reading
in the data and the current design.

Some code is shown for running the model and PF and generating
log-likelihood estimates for the design, as well as producing
some plots of the particle densities (unweighted, so only an 
approximation here) against the data for comparison, if
you wish to do this. The PF code is in the `PF.R` file, the 
truncated Skellam sampler is in the `trSkellam.R` file.
