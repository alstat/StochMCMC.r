*********
StochMCMC.r
*********
*A R package for Stochastic Gradient Markov Chain Monte Carlo*

This package is part of my master's thesis entitled **Bayesian Autoregressive Distributed Lag** *via* **Stochastic Gradient Hamiltonian Monte Carlo** or BADL-SGHMC. However as the title says, this package aims to accomodate other Stochastic Gradient MCMCs in the near future. At the latest, the following are the MCMC algorithms available:

1. Metropolis-Hasting
2. Hamiltonian Monte Carlo
3. Stochastic Gradient Hamiltonian Monte Carlo

Installation
==================
To install the package, simply run the following codes

.. code-block:: R

	library(devtools)
	install_github("alstat/StochMCMC.r")

And to load the package:

.. code-block:: R
	
	library(StochMCMC)