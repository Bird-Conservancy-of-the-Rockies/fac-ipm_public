# Simulation Study

Files and folders in this sub-directory support the simulation exercise that explores the utility of the FAC-IPM under varying types and amount of data availability. Example MCMC output is also provided, to obviate the need to actually run the simulations.

## Requirements:
* program R (https://cran.r-project.org/)
* these R packages: abind, coda, mvtnorm, parallel, RColorBrewer, rjags, vioplot, xtable
* JAGS (https://mcmc-jags.sourceforge.io/)
* GNU parallel (https://www.gnu.org/software/parallel/)

## OS notes
The code to actually run the simulations is designed for use on a Linux machine.

It may work ok on Mac OSX.

It is not set up to run on Windows OS. There are two details that will hinder the code's functionality on Windows:

* GNU parallel

    This is apparently installable on Windows, but the author has not tried this.

* the `mclapply()` function in R package `parallel`, which does not work on Windows

    This can be substituted with the `parLapply()` function, which will take the same arguments, but requires a cluster of cores to be initialized before it is called.

**However, example MCMC output has been provided with this repository, so there's really no need to run the simulations yourself.  Rather, skip #1 below and go to #2 to examine the provided output.**

## To run
1. To run the simulation study, one must make successive calls to GNU parallel.  These are shown in the file `fac-ipm_public/sim/scripts/bash/example_call_to_parallel.sh`
2. To examine MCMC output provided in this repository (or the output of your own runs from #1), gather the results and produce summaries with `fac-ipm_public/sim/scripts/R/gather_simFits_mkII.r`.  Note that this script is intended to be run interactively.

