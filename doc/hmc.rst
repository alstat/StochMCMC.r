Hamiltonian Monte Carlo
===================

Implementation of the Hamiltonian Monte Carlo sampler for Bayesian inference.

.. function:: HMC(U, K, dU, dK, init_est, d)

    Construct a ``Sampler`` object for Hamiltonian Monte Carlo sampling.

    **Arguments**

        * ``U`` : the potential energy or the negative log posterior of the parameter of interest.
        * ``K`` : the kinetic energy or the negative exponential term of the log auxiliary distribution.
        * ``dU`` : the gradient or first derivative of the potential energy ``U``.
        * ``dK`` : the gradient or first derivative of the kinetic energy ``K``.
        * ``init_est`` : the initial/starting value for the markov chain.
        * ``d`` : the dimension of the posterior distribution.

    **Value**

        Returns a ``HMC`` type object.

    **Example**
      
	In order to illustrate the modeling, the data is simulated from a simple linear regression expectation function. That is the model is given by

	.. code-block:: txt 

		y_i = w_0 + w_1 x_i + e_i,   e_i ~ N(0, 1 / a)

	To do so, let :code:`B = [w_0, w_1]' = [.2, -.9]', a = 1 / 5`. Generate 200 hypothetical data:

	.. code-block:: txt

		library(gridExtra)
		library(lattice)
		library(StochMCMC)

		set.seed(123)

		# Define data parameters
		w0 <- -.3; w1 <- -.5; stdev <- 5.; a <- 1 / stdev

		# Generate Hypothetical Data
		n <- 200;
		x <- runif(n, -1, 1);
		A <- cbind(1, x);
		B <- rbind(w0, w1);
		f <- A %*% B;
		y <- f + rnorm(n, 0, a);

		my_df = data.frame(Independent = round(x, 4), Dependent = round(y, 4));

	To view the head of the data, run the following:

	.. code-block:: txt

		head(my_df)
		#   Independent Dependent
		# 1     -0.4248   -0.2297
		# 2      0.5766   -0.5369
		# 3     -0.1820   -0.2583
		# 4      0.7660   -0.7525
		# 5      0.8809   -0.9308
		# 6     -0.9089    0.1454

	Next is to plot this data which can be done as follows:

	.. code-block:: julia

		xyplot(Dependent ~ Independent, data = my_df, type = c("p", "g"), col = "black")

	.. image:: figures/plot1.png
		:width: 80%
		:align: center
		:alt: alternate text

	In order to proceed with the Bayesian inference, the parameters of the model is considered to be random modeled by a 
	standard Gaussian distribution. That is, :code:`B ~ N(0, I)`, where :code:`0` is the zero vector. The likelihood of the data is given by,

	.. code-block:: txt

		L(w|[x, y], b) = ∏_{i=1}^n N([x_i, y_i]|w, b)

	Thus the posterior is given by,

	.. code-block:: txt

		P(w|[x, y]) ∝ P(w)L(w|[x, y], b)

	To start programming, define the probabilities

	.. code-block:: R

		# The log prior function is given by the following codes:
		logprior <- function(theta, mu = zero_vec, s = eye_mat) {
		    w0_prior <- dnorm(theta[1], mu[1], s[1, 1], log = TRUE)
		    w1_prior <- dnorm(theta[2], mu[2], s[2, 2], log = TRUE)
		    w_prior <- c(w0_prior, w1_prior)

		    w_prior %>% sum %>% return
		}

		# The log likelihood function is given by the following codes:
		loglike <- function(theta, alpha = a) {
		    yhat <- theta[1] + theta[2] * x

		    likhood <- numeric()
		    for (i in 1:length(yhat)) {
			likhood[i] <- dnorm(y[i], yhat[i], alpha, log = TRUE)
		    }

		    likhood %>% sum %>% return
		}

		# The log posterior function is given by the following codes:
		logpost <- function(theta) {
		    loglike(theta, alpha = a) + logprior(theta, mu = zero_vec, s = eye_mat)
		}

	To start the estimation, define the necessary parameters for the Metropolis-Hasting algorithm

	.. code-block:: R

		# Hyperparameters
		zero_vec <- c(0, 0)
		eye_mat <- diag(2)

	Run the MCMC:

	.. code-block:: R

		set.seed(123);
		mh_object <- MH(logpost, init_est = c(0, 0))
		chain1 <- mcmc(mh_object, r = 10000)

	Extract the estimate

	.. code-block:: R

		burn_in <- 100;
		thinning <- 10;

		# Expetation of the Posterior
		est1 <- colMeans(chain1[seq((burn_in + 1), nrow(chain1), by = thinning), ])
		est1
		# [1] -0.2984246 -0.4964463

	Setup the necessary paramters including the gradients. The potential energy is the negative logposterior given by \ 
	:code:`U`, the gradient is :code:`dU`; the kinetic energy is the standard Gaussian function given by :code:`K`, with gradient :code:`dK`. Thus,

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


