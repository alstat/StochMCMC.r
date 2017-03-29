#' Hamiltonian Monte Carlo
#' @export
#' 
#' @description 
#' This is HMC
setGeneric("mcmc", function(object, ...) {
  standardGeneric("mcmc")
})

setMethod("mcmc", signature(object = "HMC"), function(object,
	leapfrog_params = c(eps = .05, tau = 20), 
	set_seed = 123, 
	r = 1e+3) {
	
	U <- object@U
	K <- object@K
	dU <- object@dU
	dK <- object@dK
	w <- object@init_est
	d <- object@d

	eps <- leapfrog_params["eps"]
	tau <- leapfrog_params["tau"]
	
	if (w %>% is.null) {
		w <- matrix(0, d, 1)
	} else {
		w <- w
	}
	H <- function(x, p) U(x) + K(p)

	if (!(set_seed %>% is.null)) 
		set.seed(set_seed)

	out <- matrix(0, r, d)
    out[1, ] <- w
	for (i in 1:(r - 1)) {		
		w <- matrix(out[i, ])
		p <- rnorm(length(w))

		oldE <- H(w, p)
		for (j in 1:tau) {
			p <- p - (eps / 2) * dU(w)
			w <- w + eps * dK(p)
			p <- p - (eps / 2) * dU(w)
		}

		newE <- H(w, p)
		dE <- newE - oldE

		if (dE < 0) {
			out[i + 1, ] <- w
		} else if (runif(1) < exp(-dE)) {
			out[i + 1, ] <- w
		} else {
			out[i + 1, ] <- out[i, ]
		}
	}
	out %>% return
}
)
