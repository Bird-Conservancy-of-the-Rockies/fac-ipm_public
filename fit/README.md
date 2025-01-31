# Model application

Files and folders in this sub-directory support fitting the FAC-IPM to Baird's Sparrow data.

## Requirements:
* program R (https://cran.r-project.org/)
* these R packages: parallel, rjags
* JAGS (https://mcmc-jags.sourceforge.io/)

## OS notes
The code is designed for use on a Linux machine.

It should work ok on Mac OSX.

It is not set up to run on Windows OS. The main hindrance to the code's functionality on Windows is the `mclapply()` function in R package `parallel` (used in `scripts/run_mod.r`), which does not work on Windows. This can be substituted with the `parLapply()` function, which will take the same arguments, but requires a cluster of cores to be initialized before it is called.

## To run

To fit the model to the provided dataset, use Rscript to call the 'scripts/run_mod.r' script, providing the required arguments.  See the initial comments in that script, for the canonical example.

## Important note about the data

Geographic coordinates are used in the breeding season, as covariates on the zero-inflation parameter, psi.  Much of the survey data gathered by the Integrated Monitoring for Bird Conservation Regions (IMBCR) program comes from private lands.  To protect landowners, the coordinates have all been altered to mask the true locations, by rounding all locations to the nearest whole degree.

