#' Hamiltonian Monte Carlo Class
#' @export
#' 
#' @description 
#' This is HMC
setClass("HMC", 
	representation(
		U = "function",
		K = "function",
		dU = "function",
		dK = "function",
		init_est = "array",
		d = "numeric"
	)
)

HMC <- function(U, K, dU, dK, init_est, d) {
	new("HMC", U = U, K = K, dU = dU, dK = dK, init_est = init_est, d = d) 
}