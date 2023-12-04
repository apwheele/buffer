# buffer

Analysis for existence of Buffer Zone in R. Joint work with Kim Rossmo. Paper to be posted.

The script `simulation_gamma.R` has the simulation results estimating power/precision to determine whether the gamma shape paramter is greater than 1.

The script `jtc_analysis.R` has the analysis of the two serial offenders, plus function to estimate non-parametric KDE (taking into account that the density cannot be smoothed below 0) and bootstrapped point-wise confidence intervals.

![](/outputs/JTC_01V2.png)

R environment is pretty tame, just need `ggplot2` and `fitdistrplus`.

