#' Metropolis-Hasting
#' @export
#' 
#' @description 
#' This is MH
setGeneric("mcmc", function(object, ...) {
  standardGeneric("mcmc")
})

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

