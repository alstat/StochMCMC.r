#' Metropolis-Hasting Class
#' @export
#' 
#' @description 
#' This is MH
MH <- function(logposterior, proposal = default_proposal, init_est = matrix(0, 2, 1), d = 2) {
  new("MH", logposterior = logposterior, proposal = proposal, init_est = init_est, d = d)
}

sigmas <- c(1, 1)
default_proposal <- function (theta) {
  random <- numeric()
  
  for (i in 1:length(theta)) {
    random[i] <- rnorm(1, theta[i], sigmas[i])
  }
  
  random %>% return
}
