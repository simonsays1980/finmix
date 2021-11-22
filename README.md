[![R-CMD-check](https://github.com/simonsays1980/finmix/actions/workflows/r-check-package.yml/badge.svg)](https://github.com/simonsays1980/finmix/actions/workflows/r-check-package.yml) [![GPLv3 license](https://img.shields.io/badge/License-GPLv3-blue.svg)](http://perso.crans.org/besson/LICENSE.html) [![GitHub release](https://img.shields.io/github/release/Naereen/StrapDown.js.svg)](https://GitHub.com/Naereen/StrapDown.js/releases/)
[![Linux](https://svgshare.com/i/Zhy.svg)](https://svgshare.com/i/Zhy.svg) [![macOS](https://svgshare.com/i/ZjP.svg)](https://svgshare.com/i/ZjP.svg) [![Windows](https://svgshare.com/i/ZhY.svg)](https://svgshare.com/i/ZhY.svg)





# finmix
**An R package for Bayesian estimation of finite mixture distributions** 

The package uses heavily C++ code to enable high performance MCMC sampling. 
To simplify handling for users each distribution comes along with some support 
functions that create required objects and starting parameters. As a result a 
user can perform Bayesian parameter estimation in only a few lines. The following 
mixtures are available: 
* Binomial, 
* Exponential, 
* Normal, 
* Multivariate Normal, 
* Poisson,
* Condiitonal Poisson
* Student-t, and 
* Multivariate Student-t.

## Literature and implementations
The methods used in this package are based on the major literature on the Bayesian estimation 
of finite mixture distributions, namely 

*Frühwirth-Schnatter, Sylvia (2006), "Finite Mixture and Markov Switching Models", 
Springer Series in Statistics ([link](https://link.springer.com/book/10.1007/978-0-387-35768-3))*.

The code of this package is related to the `bayesf` package in `Matlab` written by Sylvia 
Frühwirth-Schnatter herself. Due to the C++ extensions in this package using `Rcpp` the 
`R` implementation is almost 200x times faster than the `Matlab` version. A single MCMC run 
with 11000 iterations is usually performed within 2-3 seconds. 

## Installation
The package can be installed directly from GitHub by using the function `install_github()` 
in the `devtools` package. The package passed all checks from R CMD check on all major 
platforms and hence, should be installable on MacOS, Windows, and Linux. Be sure that you 
installed appropriate developer tools for your platform as a C++ compiler for the source 
code is needed. 

### MacOS
For MacOS the XCode Command Line Tools are needed. You should have installed these when 
installing `R`. See the [MacOSX-FAQ](https://cran.r-project.org/bin/macosx/RMacOSX-FAQ.html#Installation-of-source-packages) 
for more information on how to install source packages on MacOS.

### Windows
For Windows the [`rtools`](https://cran.r-project.org/bin/windows/Rtools/) package is needed. 
Follow the link and install this package, if you have not installed it, yet. 

## Quick start
As a quick start choose a *Poisson* mixture with two components and perform MCMC sampling: 
```
# Load the package.
library(finmix)
# Define the finite mixture model (two components). 
f_model <- model(dist="poisson", K=2, par=list(lambda=c(.3, .7)))

# Simulate data from this mixture distribution.
f_data <- simulate(f_model)

# Define the hyperparameters for MCMC sampling.
f_mcmc <- mcmc()

# Define the prior distribution for MCMC sampling.
f_prior <- priordefine(f_data, f_model)

# Set up all parameters for MCMC sampling.
(f_data ~ f_model ~ f_mcmc) %=% mcmcstart(f_data, f_model, f_mcmc)

# Perform MCMC sampling.
f_output <- mixturemcmc(f_data, f_model, f_prior, f_mcmc)
f_output
```

To estimate the parameters you simply call the function `mcmcestimate()` 
on the output from MCMC sampling: 
```
f_estimates <- mcmcestimate(f_output)
f_estimates
```
### Further functionalities
The package comes with many auxiliary functions for plotting (e.g. `plotTraces()` 
to plot MCMC traces from MCMC outputs), for subchaining (`subseq()`), for 
swapping elements (`swapElements()`), or for relabeling components (`mcmcpermute()`). 
See the documentation for further reading. 

## Some more information
This is a package worked on for years and still not fully implemented. As it is still 
maintained by a single author, please by patient with issues. 


