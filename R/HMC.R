#' Hamiltonian Monte Carlo Class
#' @export
#' 
#' @description 
#' This is HMC
HMC <- function(U, K, dU, dK, init_est, d) {
  new("HMC", U = U, K = K, dU = dU, dK = dK, init_est = init_est, d = d) 
}