# The Journey-to-Crime Buffer Zone: Measurement Issues and Methodological Challenges

Replication code for the paper, [pre-print is available here](https://www.crimrxiv.com/pub/udrmepzn/release/1).

The script `simulation_gamma.R` has the simulation results estimating power/precision to determine whether the gamma shape parameter is greater than 1.

The script `jtc_analysis.R` has the analysis of the two serial offenders, plus function to estimate non-parametric KDE (taking into account that the density cannot be smoothed below 0) and bootstrapped point-wise confidence intervals.

![](/outputs/JTC_01V2.png)

R environment is pretty tame, just need `ggplot2` and `fitdistrplus`.

