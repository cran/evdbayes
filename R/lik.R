# LIKELIHOODS

"gevlik" <-
# Computes log-likelihood of gev model at (mu,sigma,xi). Gumbel is computed for small shape.
function(par, data, trend)
{
    nas <- !is.na(data)
    data <- data[nas]
    if(!missing(trend)) trend <- trend[nas]

    gevlik2(par = par, data = data, trend = trend)
}

"gevlik2" <-
# Computes log-likelihood of gev model at (mu,sigma,xi). Gumbel is computed for small shape.
function(par, data, trend)
{
    n <- length(data)
    if(missing(trend))
      lik <- .C(C_gevlik,
          data, n, par, dns = double(1))$dns
    else
      lik <- .C(C_gevlikt,
          data, n, par, trend, dns = double(1))$dns
    lik
}

"gpdlik" <-
# Computes log-likelihood of gpd model at (loc,scale,shape). Exponential is computed for small shape.
function(par, data, trend)
{
    nas <- !is.na(data)
    data <- data[nas]
    if(!missing(trend)) trend <- trend[nas]

    gpdlik2(par = par, data = data, trend = trend)
}

"gpdlik2" <-
# Computes log-likelihood of gpd model at (loc,scale,shape). Exponential is computed for small shape.
function(par, data, trend)
{
    n <- length(data)
    if(missing(trend))
      lik <- .C(C_gpdlik,
          data, n, par, dns = double(1))$dns
    else
      lik <- .C(C_gpdlikt,
          data, n, par, trend, dns = double(1))$dns
    lik
}

"pplik" <-
# Computes log-likelihood of Poission process model at (mu,sigma,xi) with threshold u and npy observations per block. Gumbel is computed for small shape.
function(par, data, thresh, noy, trend, exact = FALSE)
{
    n <- length(data)
    thresh <- rep(thresh, length.out = n)
    nan <- !is.na(data)
    data <- data[nan]
    thresh <- thresh[nan] 
    if(!missing(trend)) trend <- trend[nan]
    
    hd <- (data > thresh)
    data <- data[hd]
    if(length(data) == 0) stop("no data above threshold")
    if(!missing(trend)) htrend <- trend[hd]

    if(!exact) {
        thn <- seq(1, length(thresh), length = length(data))
        thresh <- thresh[thn]
        if(!missing(trend)) trend <- trend[thn]
    }
      
    pplik2(par = par, data = data, thresh = thresh,
           noy = noy, trend = trend, htrend = htrend)
}

"pplik2" <-
# Computes log-likelihood of Poission process model at (mu,sigma,xi) with threshold u and npy observations per block. Gumbel is computed for small shape.
function(par, data, thresh, noy, trend, htrend)
{   
    nh <- length(data)
    n <- length(thresh)
    if(missing(trend))
      lik <- .C(C_pplik,
          data, nh, par, thresh, n, noy, dns = double(1))$dns
    else
      lik <- .C(C_pplikt,
          data, nh, par, thresh, n, noy, trend, htrend,
                dns = double(1))$dns
    lik
}

"oslik" <-
# Computes log-likelihood of gev order statistics model at (mu,sigma,xi). Gumbel is computed for small shape.
function(par, data, trend)
{
    nas <- !apply(is.na(data), 1, all)
    data <- data[nas, ,drop = FALSE]
    if(!missing(trend)) trend <- trend[nas]
    thresh <- apply(data, 1, min, na.rm = TRUE)
    rvec <- as.integer(cumsum(rowSums(!is.na(data))))
    data <- t(data)
    data <- data[!is.na(data)]
    
    oslik2(par = par, data = data, trend = trend, thresh = thresh,
           rvec = rvec)
}

"oslik2" <-
# Computes log-likelihood of gev order statistics model at (mu,sigma,xi). Gumbel is computed for small shape.
function(par, data, trend, thresh, rvec)
{
    m <- length(thresh)
    n <- length(data)
    if(missing(trend))
      lik <- .C(C_oslik,
          data, thresh, n, m, par, dns = double(1))$dns
    else
      lik <- .C(C_oslikt,
          data, thresh, n, m, rvec, par, trend, dns = double(1))$dns
    lik
}
