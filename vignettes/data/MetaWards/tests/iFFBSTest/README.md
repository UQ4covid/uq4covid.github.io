This folder contains a simulation study to test the iFFBS filter
approach to simulate seeds for MetaWards.

Firstly run `stochModel.R` to simulate data using one
of the MW design points. Outputs are placed in the `outputs`
folder (which is created if it is not present).

Then the `iFFBSFull.R` file contains a full MCMC algorithm
that estimates parameters and hidden states based on the
simulated data. This is used to assess whether the iFFBS
code is working when fitted to detailed data. This model
fits to cumulative counts for removals and deaths, as well as
the "P" class. It uses a Gaussian discrepancy error to make
the initialisation easier and to control deviations along
the time series more strongly than the Poisson.

The `iFFBSPartial.R` contains code to run a slightly simpler
version of the model which calibrates to `DH` and `DI`
only, and requires **exact** matching. This includes tests
that aim to assess how the filter performs for **fixed**
parameters but unknown states, and examines patterns over
different design points and LADs.


