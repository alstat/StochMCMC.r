#' Metropolis-Hasting
#' @export
#' 
#' @description 
#' This is MH
#' 
setGeneric("mcmc", function(object, ...) {
  standardGeneric("mcmc")
})

setClass("MH",
         representation(
           logposterior = "function",
           proposal = "function",
           init_est = "numeric",
           d = "numeric"
         )
)

setClass("HMC", 
         representation(
           U = "function",
           K = "function",
           dU = "function",
           dK = "function",
           init_est = "numeric",
           d = "numeric"
         )
)

setClass("SGHMC", 
         representation(
           dU = "function", 
           dK = "function", 
           dKSigma = "array", 
           C_mat = "array", 
           V_mat = "array", 
           init_est = "numeric", 
           d = "numeric")
)
setMethod("mcmc", signature(object = "MH"),
  function (object, r = 1e+3) {
  out <- matrix(NA, r, object@d)
  out[1, ] <- object@init_est

  for (i in 1:(r - 1)) {
    propose <- object@proposal(out[i, ])
    probab <- exp(object@logposterior(propose) - object@logposterior(out[i, ]))

    if (runif(1) < probab){
      out[i + 1, ] <- propose
    } else {
      out[i + 1, ] <- out[i, ]
    }
  }

  out %>% return
  }
)

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
  x[1, ] = w
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
  
  for (i in 1:(r - 1)) {
    p <- rnorm(d) %>% matrix
    
    for (j in 1:tau) {
      p <- p - dU(w) * eps - C %*% solve(dKSigma) %*% p + D %*% matrix(rnorm(d))
      w <- w + dK(p) * eps
    }
    
    x[i + 1, ] <- w
  }
  
  x %>% return
}
)
