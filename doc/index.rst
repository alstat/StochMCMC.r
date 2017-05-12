****
StochMCMC.jl
****

:Author: Al-Ahmadgaid B. Asaad (alasaadstat@gmail.com | http://alstatr.blogspot.com/)
:Requires: julia releases 0.4.1 or later
:Date: |today|
:License: `MIT <https://github.com/alstat/StochMCMC.jl/blob/master/LICENSE.md>`_
:Website: https://github.com/alstat/StochMCMC.jl

(Under heavy construction, Target Finish Date on Monday 12pm - Philippine Time)

A julia package for Stochastic Gradient Markov Chain Monte Carlo. The package is part of my master's thesis entitled
**Bayesian Autoregressive Distributed Lag** *via* **Stochastic Gradient Hamiltonian Monte Carlo** or **BADL-SGHMC**,
under the supervision of **Dr. Joselito C. Magadia** of School of Statistics, University of the Philippines Diliman.
This work aims to accommodate other Stochastic Gradient MCMCs in the near future.

Installation
============
To install the package, run the following

.. code-block:: julia

    Pkg.clone("https://github.com/alstat/StochMCMC.jl")

And to load the package, run

.. code-block:: julia

    using StochMCMC

Contents
~~~~~~~~~~~~~~

.. toctree::
    :maxdepth: 2

    mh.rst
    hmc.rst
    sghmc.rst


Indices
~~~~~~~~~~

* :ref:`genindex`

Hamiltonian Monte Carlo
~~~~~~~~~~~~~~~~~~~~~~
Setup the necessary paramters including the gradients. The potential energy is the negative logposterior given by :code:`U`, the gradient is :code:`dU`; the kinetic energy is the standard Gaussian function given by :code:`K`, with gradient :code:`dK`. Thus,

.. code-block:: R

	U <- function(theta) - logpost(theta)
	K <- function(p, Sigma = diag(length(p))) (t(p) %*% solve(Sigma) %*% p) / 2
	dU <- function(theta, alpha = a, b = eye_mat[1, 1]) {
	    c(
	        - alpha * sum(y - (theta[1] + theta[2] * x)),
	        - alpha * sum((y - (theta[1] + theta[2] * x)) * x)
	    ) + b * theta
	}

	dK <- function (p, Sigma = diag(length(p))) solve(Sigma) %*% p

Run the MCMC:

.. code-block:: R

	set.seed(123)
	HMC_object <- HMC(U, K, dU, dK, c(0, 0), 2)
	chain2 <- mcmc(HMC_object, leapfrog_params = c(eps = .09, tau = 20), r = 10000)

Extract the estimate

.. code-block:: R

	est2 <- colMeans(chain2[seq((burn_in + 1), nrow(chain2), by = thinning), ])
	est2
	# [1] -0.2977521 -0.5158439

Stochastic Gradient Hamiltonian Monte Carlo
~~~~~~~~~~~~~~~~~~~~~
Define the gradient noise and other parameters of the SGHMC:

.. code-block:: R

	dU_noise <- function(theta, alpha = a, b = eye_mat[1, 1]) {
	    c(
	        - alpha * sum(y - (theta[1] + theta[2] * x)),
	        - alpha * sum((y - (theta[1] + theta[2] * x)) * x)
	    ) + b * theta + matrix(rnorm(2), 2, 1)
	}

Run the MCMC:

.. code-block:: R

	set.seed(123)
	SGHMC_object <- SGHMC(dU_noise, dK, diag(2), diag(2), diag(2), init_est = c(0, 0), 2)
	chain3 <- mcmc(SGHMC_object, leapfrog_params = c(eps = .09, tau = 20), r = 10000)

Extract the estimate:

.. code-block:: R

	est3 <- colMeans(chain3[seq((burn_in + 1), nrow(chain3), by = thinning), ])
	est3
	# [1] -0.2920243 -0.4729136

Plot it

.. code-block:: R

	p0 <- xyplot(y ~ x, type = c("p", "g"), col = "black") %>%
	    update(xlab = "x", ylab = "y")

	p1 <- histogram(chain3[, 1], col = "gray50", border = "white") %>%
	    update(xlab = expression(paste("Chain Values of ", w[0]))) %>%
	    update(panel = function (x, ...) {
	        panel.grid(-1, -1)
	        panel.histogram(x, ...)
	        panel.abline(v = w0, lty = 2, col = "black", lwd = 2)
	  })

	p2 <- histogram(chain3[, 2], col = "gray50", border = "white") %>%
	    update(xlab = expression(paste("Chain Values of ", w[1]))) %>%
	    update(panel = function (x, ...) {
	        panel.grid(-1, -1)
	        panel.histogram(x, ...)
	        panel.abline(v = w1, lty = 2, col = "black", lwd = 2)
	  })

	p3 <- xyplot(chain3[, 1] ~ 1:nrow(chain3[, ]), type = c("g", "l"), col = "gray50", lwd = 1) %>%
	    update(xlab = "Iterations", ylab = expression(paste("Chain Values of ", w[0]))) %>%
	    update(panel = function (x, y, ...) {
	        panel.xyplot(x, y, ...)
	        panel.abline(h = w0, col = "black", lty = 2, lwd = 2)
	  })

	p4 <- xyplot(chain3[, 2] ~ 1:nrow(chain3[,]), type = c("g", "l"), col = "gray50", lwd = 1) %>%
	    update(xlab = "Iterations", ylab = expression(paste("Chain Values of ", w[1]))) %>%
	    update(panel = function (x, y, ...) {
	        panel.xyplot(x, y, ...)
	        panel.abline(h = w1, col = "black", lty = 2, lwd = 2)
	  })

	p5 <- xyplot(chain3[, 2] ~ chain3[, 1]) %>%
	    update(type = c("p", "g"), pch = 21, fill = 'white', col = "black") %>%
	    update(xlab = expression(paste("Chain Values of ", w[0]))) %>%
	    update(ylab = expression(paste("Chain Values of ", w[1]))) %>%
	    update(panel = function (x, y, ...) {
	        panel.xyplot(x, y, ...)
	    })

	p6 <- xyplot(y ~ x, col = "black", fill = "gray80", cex = 1.3, type = "p", pch = 21) %>% 
	    update(xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), panel = function(x, y, ...) {
	        panel.grid(h = -1, v = -1)
	        xseq <- seq(-1, 1, length.out = 100)
	        for (i in seq((burn_in + 1), nrow(chain3), by = thinning)) {
	            yhat <- chain3[i, 1] + chain3[i, 2] * xseq
	            panel.xyplot(xseq, yhat, type = "l", col = "gray50")  
	        }
	        panel.xyplot(x, y, ...)
	        panel.xyplot(xseq, est3[1] + est3[2] * xseq, type = "l", col = "black", lwd = 2)
	  })

	acf1 <- acf(chain3[seq((burn_in + 1), nrow(chain3), by = thinning), 1], plot = FALSE)
	acf2 <- acf(chain3[seq((burn_in + 1), nrow(chain3), by = thinning), 2], plot = FALSE)
	p7 <- xyplot(acf1$acf ~ acf1$lag, type = c("h", "g"), lwd = 2, col = "black") %>%
	    update(xlab = "Lags", ylab = expression(paste("Autocorrelations of ", w[1])))

	p8 <- xyplot(acf2$acf ~ acf2$lag, type = c("h", "g"), lwd = 2, col = "black") %>%
	    update(xlab = "Lags", ylab = expression(paste("Autocorrelations of ", w[1])))

	grid.arrange(p0, p1, p2, p3, p4, p5, p6, p7, p8, ncol = 3)

.. image:: figures/plot2.png
	:width: 100%
	:align: center
	:alt: alternate text
