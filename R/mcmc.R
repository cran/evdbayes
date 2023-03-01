
## POSTERIOR DENSITIES
"dpost" <-
  ## Computes log-posterior density at (mu,sigma,xi). Includes optional normal trend for location.
  function(par, prior, lh, mix = FALSE, pMassProb, normPi0, pMass, ...)
  ## mix : logical. Is there a point mass for the prior distribution
{
  if ( missing(pMass) && mix )
    stop("pMass should be present.")
  if ( missing(pMassProb) && mix )
    stop("pMassProb should be present.")
  if ( missing(normPi0) && mix )
    stop("normPi0 should be present.")
  
  dprior <- do.call(prior$prior, c(list(par = par), prior[-1]))

  if (mix == TRUE){
    if (par[3] == pMass)
      dprior <- dprior + log(pMassProb) - log(normPi0)
    else
      dprior <- dprior + log(1- pMassProb)
  }
  
  switch(lh,
         gev = dprior + gevlik2(par, ...),
         gpd = dprior + gpdlik2(par, ...),
         pp = dprior + pplik2(par, ...),
         os = dprior + oslik2(par, ...),
         none = dprior)
}

## SAMPLER
"gibbs" <-
  function(n, init, prior, lh, ..., psd, thin, burn)
{
  if(prior$trendsd) np <- 4
  else np <- 3
  param <- c("mu","sigma","xi")
  if(prior$trendsd) param <- c(param, "mutrend")
  dpst <- function(a) dpost(a, prior, lh, ...)
  
  mc <- .Call(C_gibbs, as.integer(n), as.integer(np),
              as.integer(thin), as.double(init), as.double(psd),
              quote(dpst(x)), new.env())
  
  naccept <- mc[[2]]
  nexternal <- mc[[3]]
  mc <- matrix(mc[[1]], ncol = np, byrow = TRUE)
  dimnames(mc) <- list(seq(0, n, thin), param)
  mc <- mc[as.numeric(rownames(mc)) >= burn, , drop = FALSE]
  nexternal[np+1] <- sum(nexternal)/np
  naccept[np+1] <- sum(naccept)/np
  ar <- round(c(naccept/n, nexternal/n),2)
  ardn <- list(c("acc.rates","ext.rates"), c(param,"total"))
  ar <- matrix(ar, ncol = np+1, byrow = TRUE, dimnames = ardn)
  structure(mc, ar = ar)
}

## WRAPPER TO SAMPLER
"posterior" <-
  function(n, init, prior, lh = c("none","gev","gpd","pp","os"), ..., psd,
           burn = 0, thin = 1)
{
  if (!inherits(prior, "evprior")) 
    stop("`prior' must be a prior distribution")
  if(!prior$trendsd && (length(init) != 3 || mode(init) != "numeric"))
    stop("`init' must be a numeric vector of length three")
  if(prior$trendsd && (length(init) != 4 || mode(init) != "numeric"))
    stop("`init' must be a numeric vector of length four")
  if(init[2] <= 0)
    stop("initial value of scale parameter must be positive")
  if(!prior$trendsd && (length(psd) != 3 || mode(psd) != "numeric"))
    stop("`psd' must be a numeric vector of length three")
  if(prior$trendsd && (length(psd) != 4 || mode(psd) != "numeric"))
    stop("`psd' must be a numeric vector of length four")
  if(any(psd <= 0))
    stop("`psd'  must contain positive values")
  if(any(psd <= 0))
    stop("`psd'  must contain positive values")
  if(burn > n)
    stop("burn-in period is larger than run length")
  if(thin > n)
    stop("thinning interval is larger than run length")
  
  lh <- match.arg(lh)
  ar <- list(...)
  
  if(lh == "pp") {
    nn <- length(ar$data)
    ar$thresh <- rep(ar$thresh, length.out = nn)
    nan <- !is.na(ar$data)
    ar$data <- ar$data[nan]
    ar$thresh <- ar$thresh[nan]
    if(!is.null(ar$trend)) ar$trend <- ar$trend[nan]
    
    hd <- (ar$data > ar$thresh)
    ar$data <- ar$data[hd]
    if(length(ar$data) == 0) stop("no data above threshold")
    if(!is.null(ar$trend)) ar$htrend <- ar$trend[hd]

    if(is.null(ar$exact)) ar$exact <- FALSE
    if(!ar$exact) {
      thn <- seq(1, length(ar$thresh), length = length(ar$data))
      ar$thresh <- ar$thresh[thn]
      if(!is.null(ar$trend)) ar$trend <- ar$trend[thn]
    }
    ar$exact <- NULL
  }
  
  if(lh == "gev") {
    nas <- !is.na(ar$data)
    ar$data <- ar$data[nas]  
    if(!is.null(ar$trend)) ar$trend <- ar$trend[nas]
  }
  if(lh == "gpd") {
    nas <- !is.na(ar$data)
    ar$data <- ar$data[nas]  
    if(!is.null(ar$trend)) ar$trend <- ar$trend[nas]
  }
  if(lh == "os") {
    nas <- !apply(is.na(ar$data), 1, all)
    ar$data <- ar$data[nas, ,drop = FALSE]  
    if(!is.null(ar$trend)) ar$trend <- ar$trend[nas]
    ar$thresh <- apply(ar$data, 1, min, na.rm = TRUE)
    ar$rvec <- as.integer(cumsum(rowSums(!is.na(ar$data))))
    ar$data <- t(ar$data)
    ar$data <- ar$data[!is.na(ar$data)]
  }

  initar <- c(list(par = init, prior = prior, lh = lh), ar)
  inittest <- do.call("dpost", initar) 
  if(is.infinite(inittest))
    stop("density is zero at initial parameter values")

  ar <- c(list(n = n, init = init, prior = prior, lh = lh), ar,
          list(psd = psd, thin = thin, burn = burn))
  do.call("gibbs", ar)
}



