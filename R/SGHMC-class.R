#' Stochastic Gradient Hamiltonian Monte Carlo Class
#' @export
#' 
#' @description 
#' This is SGHMC
setClass("SGHMC", 
	representation(
		dU = "function", 
		dK = "function", 
		dKSigma = "array", 
		C_mat = "array", 
		V_mat = "array", 
		init_est = "array", 
		d = "numeric")
	)

SGHMC <- function(dU, dK, dKSigma, C_mat, V_mat, init_est, d) {
	new("SGHMC", dU = dU, dK = dK, dKSigma = dKSigma, C_mat = C_mat, V_mat = V_mat, init_est = init_est, d = d)
}