#' Stochastic Gradient Hamiltonian Monte Carlo Class
#' @export
#' 
#' @description 
#' This is SGHMC
SGHMC <- function(dU, dK, dKSigma, C_mat, V_mat, init_est, d) {
  new("SGHMC", dU = dU, dK = dK, dKSigma = dKSigma, C_mat = C_mat, V_mat = V_mat, init_est = init_est, d = d)
}