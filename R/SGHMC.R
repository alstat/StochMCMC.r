#' Stochastic Gradient Hamiltonian Monte Carlo
#' @export
#' 
#' @description 
#' This is SGHMC
setGeneric("mcmc", function(object, ...) {
  standardGeneric("mcmc")
})

setMethod("mcmc", signature(object = "SGHMC"), function(object, 
	leapfrog_params = c(eps = .05, tau = 20), 
	set_seed = 123, 
	r = 1e+3) {
	
	dU <- object@dU
	dK <- object@dK
	dKSigma <- object@dKSigma
	C <- object@C_mat
	V <- object@V_mat
	w <- object@init_est
	d <- object@d

	eps <- leapfrog_params["eps"]
	tau <- leapfrog_params["tau"]

	if (w %>% is.null) {
		w = matrix(0, d, 1)
	} else {
		w = w
	}

	x <- matrix(0, r, d)
	B <- .5 * V * eps
	D <- sqrt(2 * (C - B) * eps)

	if ((dim(B)[1] != dim(C)[1]) & (dim(B)[2] != dim(C)[2])) {
		"C and V should have the same dimension." %>% error
	} else {
		if (sum(dim(B)) > 1) {
			if (det(B) > det(C)) {
				"eps is too big. Consider decreasing it." %>% error	
			}
		} else {
			if (B > C) {
				"eps is too big. Consider decreasing it." %>% error
			}
		}
	}

	for (i in 1:r) {
		p <- rnorm(d) %>% matrix

		for (j in 1:tau) {
			p <- p - dU(w) * eps - C %*% solve(dKSigma) %*% p + D %*% matrix(rnorm(d))
			w <- w + dK(p) * eps
		}

		x[i, ] <- w
	}

	x %>% return
}
)