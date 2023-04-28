# ===========================================================================
# File: "main.R"
#                        Created: 2017-05-30 09:18:44
#              Last modification: 2022-02-01 11:30:37
# Author: Bernard Desgraupes
# e-mail: <bernard.desgraupes@parisnanterre.fr>
# This file is part of the boodd project.
# ===========================================================================

## 
 # ------------------------------------------------------------------------
 # 
 # "boots(x,func,B,...)" --
 # 
 # Bootstrap for the iid case.
 # 
 # ------------------------------------------------------------------------
 ##
boots <- function(x,func,B,smooth=FALSE,moonsize=NULL,mreplace=TRUE,...) {
	# Test the value returned by func
	y <- func(x,...)
	if (!is.vector(y)) {
		stop("Function 'func' must return a vector")
	}
	len <- length(y)
	cnames <- names(y)
	res <- matrix(nrow=B,ncol=len)
	if (is.matrix(x)) {
		n <- nrow(x)
	} else {
		n <- length(x)
	}
	# Moon ('m out of n') bootstrap parameters
	if (is.null(moonsize)) {
		m <- n
		repl <- TRUE
	} else {
		m <- moonsize
		repl <- mreplace
	}
	hn <- 0
	if (smooth) {
		# Estimate unbiased cross-validation bandwith for smoothing
		s <- sd(as.vector(x))
		hn = bw.ucv(as.vector(x),lower=s*0.000001,upper=s*2)
	} 
	if (is.matrix(x)) {
		# N is different for each column
		num <- m*ncol(x)
	} else {
		num <- m
	}
	N <- rep(0,num)
   	for (i in 1:B) {
		ind <- sample(1:n,size=m,replace=repl)
		if (smooth) {
			N <- rnorm(num,0,hn)
		}
		if (is.matrix(x)) {
			res[i,] <- func(x[ind,]+matrix(N,nrow=m),...)
		} else {
			res[i,] <- func(x[ind]+N,...)
		}
	}
	if (len == 1) {
		res <- as.vector(res)
	} else if (!is.null(cnames)) {
		colnames(res) <- cnames
	}
	obj <- list(s=res,Tn=y)
	class(obj) <- "boodd"
	attr(obj,"kind") <- "iid"
	attr(obj,"func") <- func
    return(obj)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "plot.boodd <- function(x,...)" --
 # 
 # Plot an object of class 'boodd' returned by the
 # boots/bootsemi/blockboot/regenboot/etc functions. 
 # See R-3.2.2/src/library/stats/R/plot.lm.R
 # 
 # ------------------------------------------------------------------------
 ##
plot.boodd <- function(x,with.density=TRUE,which,byrow=FALSE,...) {
	iargs <- list(...)
	s <- x$s
	if (!is.matrix(x$s)) {
		s <- matrix(x$s,ncol=1)
		colnames(s) <- deparse(substitute(x))
	}
	nc <- ncol(s)
	if (missing(which)) {
		which <- 1:nc
	}
	if (length(which) > 6) {
		which <- which[1:6]
		warning("can't display more than 6 columns.")
	}
	len <- length(which)
	if (len > 1) {
		vals <- 1:len
		if ((len %% 2) == 1) {vals <- append(vals,len+1)}
		layout(matrix(vals,nrow=2,byrow=byrow))
	}
	res <- list()
	names <- colnames(s)
	for (i in 1:len) {
		oargs <- iargs
		ci <- which[i]
		if (ci > nc) {
			stop("column index out of bounds")
		}
		X <- s[,ci]
		oargs$x <- X
		oargs$prob <- TRUE
		cname <- names[ci]
		if (is.null(cname) || cname == "") {
			cname <- paste("column",ci)
		} 
		if (length(iargs) == 0 || is.null(iargs$main) || is.na(iargs$main[i])) {
			oargs$main <- paste("Histogram of",cname)
		} else {
			oargs$main <- iargs$main[i]
		}
		if (length(iargs) != 0 && !is.null(iargs$xlab) && !is.na(iargs$xlab[i])) {
			oargs$xlab <- iargs$xlab[i]
		} else if (!is.null(cname)) {
			oargs$xlab <- cname
		} else {
			oargs$xlab <- paste("x",i,sep="")
		}
				
		if (with.density) {
			H <- hist(X,plot=FALSE)
			D <- density(X)
			yl <- c(0,max(max(H$density),max(D$y)))
			oargs$ylim <- yl
			H <- do.call(hist,oargs)
			lines(density(X))
		} else {
			H <- do.call(hist,oargs)
		}
		res[[i]] <- H
	}
	if (len == 1) {res <- res[[1]]}
	layout(1)
	invisible(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "summary.boodd <- function(object,...)" --
 # 
 # Print a summary of an object of class 'boodd' returned by
 # the boots/bootsemi/blockboot/regenboot/etc functions.
 # 
 # ------------------------------------------------------------------------
 ##
summary.boodd <- function(object,...) {
	message("kind = ",paste(attr(object,"kind"),collapse=" "))
	summary.default(object$s)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "confint.boodd <- function(object,alpha,method,recenter,...)" --
 # 
 # Calculate a confidence interval for an object of class 'boodd' returned by
 # the boots/bootsemi/blockboot/regenboot/etc functions.
 # 
 # ------------------------------------------------------------------------
 ##
confint.boodd <- function(object,alpha=0.05,method=c("perc","bperc","aboot","tboot","tsymboot","all"),recenter,...) {
	if (missing(method)) {
		method <- c("perc")
	}
	method <- match.arg(method)
	kind <- attr(object,"kind")
	if (missing(recenter)) {
		if (kind[1] == "block") {
			recenter <- (kind[2] == "movingblock" || kind[2] == "nonoverlapping")
		} else {
			recenter <- FALSE
		}
	} else if ((method != "tboot") && (method != "tsymboot") && (method != "all")) {
		warning("the 'recenter' argument is used only by the 'tboot' and 'tsymboot' methods")
	}
	s <- object$s
	if (!is.matrix(s)) {
		s <- matrix(s,ncol=1)
	} 
	nc <- ncol(s)
	B <- nrow(s)
	T <- object$Tn
	wantAll <- (method == "all")
	if (wantAll) {
		methods <- c("perc","bperc","aboot","tboot","tsymboot")
		L <- list()
	} else {
		methods <- method
	}
	
	for (meth in methods) {
		if ((meth == "tboot") || (meth == "tsymboot")) {
			if ((nc %% 2) != 0) {
				msg <- paste("The",meth,"method expects an even number of columns.")
				if (wantAll) {
					warning(msg)
					res <- matrix(nrow=1,ncol=2)
				} else {
					stop(msg)
				}
			} else {
				k <- nc/2
				res <- matrix(nrow=k,ncol=2)
				for (i in 1:k) {
					b <- s[,i]
					V <- s[,i+k]
					if (meth == "tboot") {
						t <- sort((b-T[i])/sqrt(V))
						q1 <- valueForPos(t,(B+1)*alpha/2)
						q2 <- valueForPos(t,(B+1)*(1-alpha/2))
						res[i,] <- T[i]-c(q2,q1)*sqrt(T[i+k])
					} else {
						# Symmetric interval using the absolute values
						t <- sort(abs(b-T[i])/sqrt(V))
						q <- valueForPos(t,(B+1)*(1-alpha/2))
						r <- q*sqrt(T[i+k])
						res[i,] <- T[i]+c(-r,r)
					}
				}
				if (recenter) {
					C <- colMeans(s[,1:k,drop=FALSE])
					R <- rowMeans(res)
					res <- res + cbind(C-R,C-R)
				}
				rownames(res) <- colnames(s[,1:k])
			}
		} else {
			res <- matrix(nrow=nc,ncol=2)
			for (i in 1:nc) {
				if (meth == "bperc") {
					S <- sort(s[,i])
					lo <- valueForPos(S,(B+1)*alpha/2)
					up <- valueForPos(S,(B+1)*(1-alpha/2))
					res[i,] <- c(lo,up)
				} else if (meth == "perc") {
					S <- sort(s[,i])
					lo <- valueForPos(S,(B+1)*alpha/2)
					up <- valueForPos(S,(B+1)*(1-alpha/2))
					res[i,] <- c(2*T[i]-up,2*T[i]-lo)
				} else if (meth == "aboot") {
					sdboot <- sd(s[,i])
					u <- qnorm(1-alpha/2)
					res[i,] <- T[i]+u*sdboot*c(-1,1)
				} 
			}
		}
		colnames(res) <- c("lwr","upr")
		if (wantAll) {
			L[[meth]] <- res
		}
	}
	
	if (wantAll) {
		return(L)
	} else {
		return(res)
	}
}


## 
 # ------------------------------------------------------------------------
 # 
 # "bootsemi(x,func,B,...,model,params,model.fit,model.sim)" --
 # 
 # Bootstrap for time series.
 # 
 # ------------------------------------------------------------------------
 ##
bootsemi <- function(x,func,B,...,model=c("ARIMA","GARCH"),params,model.fit=NULL,model.sim=NULL) {
	# Test the value returned by func
	y <- func(x,...)
	if (!is.vector(y)) {
		stop("Function 'func' must return a vector")
	}
	len <- length(y)
	cnames <- names(y)
	
	if (!is.null(model.fit) || !is.null(model.sim)) {
		if (!is.null(model.fit) && !is.null(model.sim)) {
			res <- bootsemi.generic(x,func,B,model.fit,model.sim,params,...)
		} else {
			stop("both arguments 'model.fit' and 'model.sim' are required.")
		}
	} else {
		model <- match.arg(model)
		if (model[1] == "ARIMA") {
			res <- bootsemi.arima(x,func,B,params,...)
		} else if (model[1] == "GARCH") {
			res <- bootsemi.garch(x,func,B,params,...)
		} else {
			stop("unknown model ",model[1])
		}
	}
	
	if (len == 1) {
		res <- as.vector(res)
	} else if (!is.null(cnames)) {
		colnames(res) <- cnames
	}
	obj <- list(s=res,Tn=y)
	class(obj) <- "boodd"
	attr(obj,"kind") <- "semi"
	attr(obj,"func") <- func
    return(obj)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "bootsemi.generic <- function(x,func,B,model.fit,model.sim,params,...)" --
 # 
 # Bootstrap for generically modelled time series.
 # 
 # ------------------------------------------------------------------------
 ##
bootsemi.generic <- function(x,func,B,model.fit,model.sim,params,...) {
    n <- length(x)
	y <- func(x,...)
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	
	# Fit the model
	fit <- model.fit(x,params)
	# Centered residuals epsilon-tilde
	eps_hat <- fit$residuals
	eps_tilde <- eps_hat - mean(eps_hat)
	
	for (i in 1:B) {
		# Bootstrap the centered residuals
		eps_star <- sample(eps_tilde,size=n,replace=TRUE)
		xstar <- model.sim(fit,eps_star,params)
		res[i,] <- func(xstar,...)
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "bootsemi.arima <- function(x,func,B,params,...)" --
 # 
 # Bootstrap for arima modelled time series. The arima() function from the 
 # stats package uses the following definition of the ARMA model:
 #     X[t] = a[1]X[t-1] + … + a[p]X[t-p] + e[t] + b[1]e[t-1] + … + b[q]e[t-q]
 # 
 # ------------------------------------------------------------------------
 ##
bootsemi.arima <- function(x,func,B,params,...) {
    n <- length(x)
	y <- func(x,...)
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	lor <- length(params)
	if ((lor != 2) && (lor != 3)) {
		stop("wrong 'params' argument. Should be vector (p,q) or (p,d,q).")
	}
	if (length(params) == 2) {
		# Assume (p,q) specified and d=0
		params[3] <- params[2]
		params[2] <- 0
	}
	p <- params[1]
	d <- params[2]
	q <- params[3]
	phi <- 0
	theta <- 0
	# Fit the ARIMA model
	fit <- arima(x,params)
	# Estimated coefficients
	coeffs <- coef(fit)
	if (p > 0) {
		phi <- coeffs[1:p]
	}
	if (q > 0) {
		theta <- coeffs[(p+1):(p+q)]
	}	
	
	# Centered residuals epsilon-tilde
	eps_hat <- residuals(fit)
	eps_tilde <- eps_hat - mean(eps_hat)
	xstar <- numeric(n)
	if (p > 0) {
		xstar[1:p] <- x[1:p]
	}
	
	for (i in 1:B) {
		# Bootstrap the centered residuals
		neweps <- sample(eps_tilde,size=n,replace=TRUE)
		for (t in (p+1):n) {
			xstar[t] <- neweps[t]
			if (p>0) {
				xstar[t] <- xstar[t]+sum(phi*xstar[(t-1):(t-p)])
			}
			if (q>0 && t>q) {
				xstar[t] <- xstar[t]+sum(theta*neweps[(t-1):(t-q)])
			}
		}
		res[i,] <- func(xstar,...)
	}
	return(res)
}

# tools::assertError

## 
 # ------------------------------------------------------------------------
 # 
 # "bootsemi.garch <- function(x,func,B,params,...)" --
 # 
 # Residual bootstrap for GARCH modelled variance. See Shimizu p.69sq.
 # 
 # epsilon_t = sqrt(h_t) * eta_t = sigma_t * eta_t
 # h_t = sigma_t^2 = omega + \sum _{i=1}^{q} alpha_{i} epsilon_{t-i}^{2} + \sum_{i=1}^{p} beta_{i} sigma_{t-i}^{2}
 # 
 # ------------------------------------------------------------------------
 ##
bootsemi.garch <- function(x,func,B,params,...) {
	n <- length(x)
	y <- func(x,...)
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	lor <- length(params)
	if ((lor != 1) && (lor != 2)) {
		stop("wrong 'params' argument. Should be vector (q) or (p,q).")
	}
	if (length(params) == 1) {
		# Assume it is an ARCH(q)
		params[2] <- params[1]
		params[1] <- 0
	}
	p <- params[1]
	q <- params[2]
	if (q < 1) {stop("One must have q>0.")}
	
	# Step 1: fit the GARCH model
	fit <- garch(x,params)
	coeffs <- coef(fit)
	omega <- coeffs[1]
	alpha <- coeffs[2:(q+1)]
	beta <- coeffs[(q+2):(p+q+1)]
	if (alpha+beta >= 1) {
		warning("the fitted series does not verify the stationarity condition (the sum of the coefficients is greater than 1)")
	}
	
	# Step 2: compute the estimated heteroscedasticity. It is returned in 
	# the first column of fit$fitted.values. NB: the garch function puts 
	# NAs in the first maxpq values. We restore them to the unconditional 
	# variance (which is the value they receive in the tseries_pred_garch 
	# function. See source code of the tseries package).
		
	# Unconditional variance
	uv <- omega/(1-alpha-beta)
	sigt <- fit$fitted.values[,1]
	sigt[is.na(sigt)] <- sqrt(uv)
	
	# Step 3: the residuals are not yet standardised
	etahat <- as.vector(x)/sigt
	etatil <- scale(etahat)
	
	for (i in 1:B) {
		# Step 4: draw from eta_tilde
		etastar <- sample(etatil,n,replace=TRUE)
		
		# Step 5: generate the bootstrapped process
		epsstar <- sigt*etastar
		res[i,] <- func(epsstar,...)
	}
	
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "sieveboot <- function(x,func,B,order,...)" --
 # 
 # Autoregressive Sieve Bootstrap for stationary time series that admit a 
 # MA(∞) representation.
 # 
 # ------------------------------------------------------------------------
 ##
sieveboot <- function(x,func,B,order=NULL,...) {
	x <- as.vector(x)
    n <- length(x)
	y <- func(x,...)
	len <- length(y)
	cnames <- names(y)
	res <- matrix(nrow=B,ncol=len)
	
	if (is.null(order)) {
		#ln <- floor(n^0.25/sqrt(log(n)))
		ln <- floor(4*n^0.25/log(n)^0.5)
	} else {
		ln <- order
	}
	if (n <= ln) {
		stop("length of series must be greater than order")
	}
	
	# Fit the linear autoregressive sieve process
	fit <- arima(x, order=c(ln,0,0))
	coeffs <- coef(fit)
	phi <- coeffs[1:ln]
	eps_hat <- residuals(fit)
	# Remove the first ln values and center on the mean
	eps_hat <- eps_hat[-(1:ln)]
	eps_tilde <- eps_hat - sum(eps_hat)/(n-ln)

	xstar <- numeric(n)
	xstar[1:ln] <- x[1:ln]
	
	for (i in 1:B) {
		# Bootstrap the centered residuals
		eps_star <- sample(eps_tilde,size=n,replace=TRUE)
		for (t in (ln+1):n) {
			xstar[t] <- sum(phi*xstar[(t-1):(t-ln)])+eps_star[t]
		}
		res[i,] <- func(xstar,...)
	}
	if (len == 1) {
		res <- as.vector(res)
	} else if (!is.null(cnames)) {
		colnames(res) <- cnames
	}
	obj <- list(s=res,Tn=y)
	class(obj) <- "boodd"
	attr(obj,"kind") <- "sieve"
	attr(obj,"func") <- func
    return(obj)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "blockboot(x,func,B,length.block,method,period,...)" --
 # 
 # Block bootstrap for time series. Method names maybe abbreviated.
 # 
 # ------------------------------------------------------------------------
 ##
blockboot <- function(x,func,B,length.block=NULL,...,method=c("movingblock","nonoverlapping","circular","stationary","seasonal"),period) {
	method <- match.arg(method)
	if (is.null(length.block)) {
		defaultlen <- length(x)^(1/3)
		if (method[1] == "seasonal") {
			length.block <- ceiling(defaultlen)
		} else {
			if (method[1] == "circular") {
				# BstarCB
				blstar <- b.star(x)[2]
			} else {
				# BstarSB
				blstar <- b.star(x)[1]
			}
			length.block <- ceiling(max(c(blstar,defaultlen)))
		}
	}
	
	if (length.block >= length(x)) {
		stop("the block length must be less than the size of the time series")
	}
	# Test the value returned by func
	y <- func(x,...)
	if (!is.vector(y)) {
		stop("Function 'func' must return a vector")
	}
	len <- length(y)
	cnames <- names(y)

	if (method[1] == "movingblock") {
		res <- blockboot.mb(x,func,B,length.block,...)
	} else if (method[1] == "nonoverlapping") {
		res <- blockboot.nooverlap(x,func,B,length.block,...)
	} else if (method[1] == "circular") {
		res <- blockboot.circular(x,func,B,length.block,...)
	} else if (method[1] == "stationary") {
		res <- blockboot.stationary(x,func,B,length.block,...)
	} else if (method[1] == "seasonal") {
		if (missing(period)) {
			if (is.ts(x)) {
				warning("Missing 'period' argument, using the frequency of the series (",frequency(x),").")
				period <- frequency(x)
			} else {
				stop("missing 'period' argument")
			}
		} 
		res <- blockboot.seasonal(x,func,B,length.block,period,...)
	} 
	
	if (len == 1) {
		res <- as.vector(res)
	} else if (!is.null(cnames)) {
		colnames(res) <- cnames
	}
	obj <- list(s=res,Tn=y)
	class(obj) <- "boodd"
	attr(obj,"kind") <- c("block",method[1])
	attr(obj,"func") <- func
    return(obj)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "blockboot.mb <- function(x,func,B,length.block,...)" --
 # 
 # Block bootstrap for time series using moving blocks.
 # 
 # ------------------------------------------------------------------------
 ##
blockboot.mb <- function(x,func,B,length.block,...) {
    n <- length(x)
	y <- func(x,...)
	b <- length.block
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	
	# Number of blocks
	q <- n-b+1
	# Number of sections
	s <- ((n-1) %/% b)+1
	# Bootstrapped series
	nx <- numeric(s*b)
	for (i in 1:B) {
		# Bootstrap the blocks (represented by their starting index)
		bb <- sample(1:q,size=s,replace=TRUE)
		pos <- 1
		# Rebuild a new series with the chosen blocks
		for (idx in bb) {
			nx[pos:(pos+b-1)] <- x[idx:(idx+b-1)]
			pos <- pos+b
		}
		res[i,] <- func(nx[1:n],...)
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "blockboot.nooverlap <- function(x,func,B,length.block,...)" --
 # 
 # Block bootstrap for time series using non overlapping blocks.
 # 
 # ------------------------------------------------------------------------
 ##
blockboot.nooverlap <- function(x,func,B,length.block,...) {
    n <- length(x)
	y <- func(x,...)
	b <- length.block
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	
	# Number of blocks
	q <- n %/% b
	# Number of sections
	s <- ((n-1) %/% b)+1
	# Bootstrapped series
	nx <- numeric(s*b)
	for (i in 1:B) {
		# Bootstrap the blocks (represented by their index)
		bb <- sample(1:q,size=q+1,replace=TRUE)
		pos <- 1
		# Rebuild a new series with the chosen blocks
		for (idx in bb) {
			start <- (idx-1)*b+1
			nx[pos:(pos+b-1)] <- x[start:(start+b-1)]
			pos <- pos+b
		}
		res[i,] <- func(nx[1:n],...)
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "blockboot.circular <- function(x,func,B,length.block,...)" --
 # 
 # Block bootstrap for time series using circular blocks.
 # 
 # ------------------------------------------------------------------------
 ##
blockboot.circular <- function(x,func,B,length.block,...) {
    n <- length(x)
	y <- func(x,...)
	b <- length.block
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	# Ensure circularity
	x <- c(x,x[1:(b-1)])
	
	# Number of sections
	s <- ((n-1) %/% b)+1
	# Bootstrapped series
	nx <- numeric(s*b)
	for (i in 1:B) {
		# Bootstrap the blocks (represented by their starting index)
		bb <- sample(1:n,size=s,replace=TRUE)
		pos <- 1
		# Rebuild a new series with the chosen blocks
		for (idx in bb) {
			nx[pos:(pos+b-1)] <- x[idx:(idx+b-1)]
			pos <- pos+b
		}
		res[i,] <- func(nx[1:n],...)
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "blockboot.stationary <- function(x,func,B,length.block,...)" --
 # 
 # Block bootstrap for time series using random length blocks taken from
 # random positions. Lengths are drawn from a geometric distribution with
 # mean equal to 'length.block' (p=1/b) and positions are drawn uniformly.
 # 
 # ------------------------------------------------------------------------
 ##
blockboot.stationary <- function(x,func,B,length.block,...) {
    n <- length(x)
	y <- func(x,...)
	b <- length.block
	p <- 1/b
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	# Ensure circularity
	x <- c(x,x)
	
	for (i in 1:B) {
		# Generate block lengths with geometric distribution
		L <- numeric(0)
		tot <- 0
		while (tot < n) {
			nl <- rgeom(1,p)
			if (nl > 0) {
				L <- c(L,nl)
				tot <- tot + nl
			}
		}
		# Number of blocks
		K <- length(L)
		# Generate block starts with uniform distribution
		S <- sample(1:n,K,replace=TRUE)
		
		# Bootstrapped series
		nx <- numeric(tot)
		pos <- 1
		# Rebuild a new series with the chosen blocks
		for (j in 1:K) {
			start <- S[j]
			b <- L[j]
			nx[pos:(pos+b-1)] <- x[start:(start+b-1)]
			pos <- pos+b
		}
		res[i,] <- func(nx[1:n],...)
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "blockboot.seasonal <- function(x,func,B,length.block,period,...)" --
 # 
 # Generalized seasonal block bootstrap (GSBB) for time series.
 # 
 # ------------------------------------------------------------------------
 ##
blockboot.seasonal <- function(x,func,B,length.block,period,...) {
	n <- length(x)
	y <- func(x,...)
	b <- length.block
	# Trim the series to an entire number of periods
	w <- n%/%period
	x <- x[1:(w*period)]
	# Ensure circularity
	x <- c(x,x[1:b])
	# Allocate the result object
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	# Number of sections
	s <- ((n-1) %/% b)+1
	# Bootstrapped series
	nx <- numeric(s*b)
	
	for (i in 1:B) {
		c <- 1
		pos <- 1
		# The number of periodic indices c, c+d, c+2*d, ... up to n is w
		base <- (0:(w-1))*period
		while (pos <= n) {
			# Draw an index
			idx <- sample(c+base,1)
			nx[pos:(pos+b-1)] <- x[idx:(idx+b-1)]
			pos <- pos+b
			c <- (c+b-1)%%period+1
		}
		res[i,] <- func(nx[1:n],...)
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "regenboot(x,func,B,...,atom,s=mean(x),eps)" --
 # 
 # Regenerative bootstrap for Markov chains.
 # 
 # ------------------------------------------------------------------------
 ##
regenboot <- function(x,func,B,...,atom,small=NULL,s=median(x),plotIt=FALSE) {
	# Test the value returned by func
	y <- func(x,...)
	if (!is.vector(y)) {
		stop("Function 'func' must return a vector")
	}
	len <- length(y)
	cnames <- names(y)
	
	if (!missing(atom)) {
		# This is the atomic case
		res <- regenboot.atom(x,func,B,atom=atom,...)
	} else {
		if (is.null(small)) {
			small <- findBestEpsilon(x,s)
		}
		res <- regenboot.smallEnsemble(x,func,B,small,...,plotIt=plotIt)
	}
	
	if (len == 1) {
		res <- as.vector(res)
	} else if (!is.null(cnames)) {
		colnames(res) <- cnames
	}
	obj <- list(s=res,Tn=y)
	class(obj) <- "boodd"
	attr(obj,"kind") <- "regenerative"
	attr(obj,"func") <- func
    return(obj)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "regenboot.atom <- function(x,func,B,atom,...)" --
 # 
 # Atomic regenerative bootstrap for finite states Markov chains.
 # 
 # ------------------------------------------------------------------------
 ##
regenboot.atom <- function(x,func,B,atom,...) {
	if (length(atom) != 1) {
		stop("atom must be a single value")
	}
	if ( !(atom %in% x) ) {
		stop("atom must be an element of x")
	}
    n <- length(x)
	y <- func(x,...)
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)

	idx <- which(x==atom)
	li <- length(idx)
	starts <- idx+1
	lens <- diff(c(idx,n))
	# Drop the last block, either because it is incomplete or because it is
	# past the last position
	nb <- li - 1
	
	for (i in 1:B) {
		nx <- numeric(n+max(lens)-1)
		tot <- 0
		pos <- 1
		while (pos <= n) {
			# Draw a block
			bi <- sample(1:nb,1)
			start <- starts[bi]
			len <- lens[bi]
			nx[pos:(pos+len-1)] <- x[start:(start+len-1)]
			pos <- pos + len
		}
		res[i,] <- func(nx[1:n],...)
	}
		
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "regenboot.smallEnsemble <- function(x,func,B,small,...)" --
 # 
 # Small ensemble regenerative bootstrap for homogeneous Markov chains.
 # 
 # ------------------------------------------------------------------------
 ##
regenboot.smallEnsemble <- function(x,func,B,small,...,plotIt=FALSE) {
	if (class(small)[1] != "smallEnsemble") {
		stop("expecting an object of class 'smallEnsemble'.")
	}
	
	y <- func(x,...)
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)

	s <- small$s
	eps <- small$epsilon
	delta <- small$delta
	p_XiXip1 <- small$trans	
	n <- length(x)

	# Find indices of pairs (X_i,X_{i+1}) belonging to the small ensemble
	isInSmall <- (x[1:(n-1)]>=s-eps)*(x[1:(n-1)]<=s+eps)*(x[2:n]>=s-eps)*(x[2:n]<=s+eps)

	# Compute the Bernoulli parameter if (X_i,X_{i+1}) \in S
	prob_regen <- delta*isInSmall/p_XiXip1

	# Simulate the Bernoulli drawing
	regen <- c((prob_regen>runif(n-1)),0)

	# Build the partition in blocks (ignoring first and last blocks)
	blocknums <- cumsum(c(0,regen[1:(n-1)]))
	if (regen[n-1]==1) {
		nb <- max(blocknums)
	} else {
		nb <- max(blocknums)-1
	}
	data <- cbind(1:n,blocknums)
	# S: start index of block - L: block length
	S <- numeric(nb)
	L <- numeric(nb)
	for (i in 1:nb) {
		aux <- subset(data,data[,2]==i,1)
		S[i] <- aux[1]
		L[i] <- nrow(aux)
	}
	
	# Bootstrap now
	for (i in 1:B) {
		nx <- numeric(n+max(L)-1)
		tot <- 0
		pos <- 1
		while (pos <= n) {
			# Draw a block
			bi <- sample(1:nb,1)
			start <- S[bi]
			len <- L[bi]
			nx[pos:(pos+len-1)] <- x[start:(start+len-1)]
			pos <- pos + len
		}
		res[i,] <- func(nx[1:(n-1)],...)
	}

	if (plotIt) {
		# Display the series
		plot(x,type="l",col="gray",xlab="",ylab="x",pch=16,main="(a)RBB Blocks")
		abline(h=s,col="gray")
		abline(h=c(s-eps,s+eps),lty=2,col="red")

		#  Draw the blocks
		for (i in 1:nb) {
			b <- S[i]
			l <- L[i]
			lines(b:(b+l),x[b:(b+l)],col=i)
			abline(v=b,lty=3,col="blue")
		}
	}
		
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "freqboot(x,XI,g,B,kernel="normal")" --
 # 
 # Bootstrap in frequency domain (FDB).
 # 
 # ------------------------------------------------------------------------
 ##
freqboot <- function(x,XI,g,B,kernel="normal",bandwidth) {
	n <- length(x)
	p <- length(XI)
	# Check the arguments
	for (i in 1:p) {
		if (!is.list(XI) | !is.function(XI[[i]])) {
			stop("XI must be a list of functions")
		}
	}
	if (length(formalArgs(g)) != p) {
		stop("g must be a numeric function with as many arguments as the length of XI (",p,")")
	}
	y <- do.call(g,as.list(rep(1,p)))
	if (!is.vector(y)) {
		stop("Function 'g' must return a vector")
	}
# 	if (missing(bandwidth)) {h <- bandw1(x)}
	if (missing(bandwidth)) {h <- sd(x)*n^(-1/3)}

	# Periodograms
	x <- ts(x,frequency=1)
	# Compute via FFT
	P <- spec.pgram(x,plot=FALSE,taper=0,fast=FALSE,detrend=FALSE)
	# The returned object P is scaled with factor 1/frequency(x) and without
	# the 1/(2*pi) prefix. It contains the periodograms for the positive
	# Fourier frequencies
	specs <- P$spec/(2*pi)
	freqs <- P$freq*2*pi
	n0 <- length(specs)
	
	# Precompute the \xi_i(\lambda_{jn})
	xiljn <- matrix(nrow=p,ncol=n0)
	for (i in 1:p) {
		xiljn[i,] <- XI[[i]](freqs)
	}

	# Precompute the smoothing coefficients
	smc <- smoothingCoefficients(n,h,kernel)

	# Initial statistic
	V <- numeric(p)
	for (i in 1:p) {
		V[i] <- mean(xiljn[i,]*specs)
	}
	y <- do.call(g,as.list(V))
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	cnames <- names(y)

	# Extend periodograms and frequencies to ]-pi,pi[ 
	P0 <- sum(x)^2/(2*pi*n)
	In <- c(rev(specs),P0,specs)
	ljn <- c(-rev(freqs),0,freqs)

	# Compute \hat{f_n} for positive Fourier frequencies 
	fnhat <- numeric(n0)
	for (i in 1:n0) {
		K <- smc[(n0-i+1):(3*n0-i+1)]
		fnhat[i] <- 2*pi*mean(K*In)/h
	}

	# Compute \hat{\epsilon_jn}
	eps_hat <- specs/fnhat
	# Normalize by the mean
	eps_tilde <- eps_hat/mean(eps_hat)
	# Bootstrap
	for (j in 1:B) {
		# Draw epsilon values randomly with replacement
		eps_star <- sample(eps_tilde,n0,replace=TRUE)
		I_star <-fnhat*eps_star	
		for (i in 1:p) {
			V[i] <- mean(xiljn[i,]*I_star)
		}
		res[j,] <- do.call(g,as.list(V))
	}

	if (len == 1) {
		res <- as.vector(res)
	} else if (!is.null(cnames)) {
		colnames(res) <- cnames
	}
	obj <- list(s=res,Tn=y)
	class(obj) <- "boodd"
	attr(obj,"kind") <- "frequency"
	attr(obj,"func") <- g
    return(obj)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "aidedboot(x,XI,g,B,order,kernel,bandwidth)" --
 # 
 # Aided Frequency Bootstrap (AFB).
 # 
 # ------------------------------------------------------------------------
 ##
aidedboot <- function(x,XI,g,B,order=NULL,kernel="normal",bandwidth) {
	n <- length(x)
	p <- length(XI)
	# Check the arguments
	for (i in 1:p) {
		if (!is.list(XI) | !is.function(XI[[i]])) {
			stop("XI must be a list of functions")
		}
	}
	if (length(formalArgs(g)) != p) {
		stop("g must be a numeric function with as many arguments as the length of XI (",p,")")
	}
	y <- do.call(g,as.list(rep(1,p)))
	if (!is.vector(y)) {
		stop("Function 'g' must return a vector")
	}
	if (missing(bandwidth)) {h <- bandw1(x)}
	if (is.null(order)) {
		#ln <- floor(n^0.25/sqrt(log(n)))
		ln <- floor(4*(n/log(n))^0.25)
	} else {
		ln <- order
	}
	if (n <= ln) {
		stop("length of series must be greater than order")
	}

	# Initial periodograms
	x <- ts(x,frequency=1)
	P <- spec.pgram(x,plot=FALSE,taper=0,fast=FALSE,detrend=FALSE)
	specs <- P$spec/(2*pi)
	freqs <- P$freq*2*pi
	n0 <- length(specs)
	P0 <- sum(x)^2/(2*pi*n)
	I_n <- c(rev(specs),P0,specs)
	ljn <- c(-rev(freqs),0,freqs)

	# Precompute the \xi_i(\lambda_{jn})
	xiljn <- matrix(nrow=p,ncol=n0)
	for (i in 1:p) {
		xiljn[i,] <- XI[[i]](freqs)
	}

	# Precompute the smoothing coefficients
	smc <- smoothingCoefficients(n,h,kernel)

	# Initial statistic
	V <- numeric(p)
	for (i in 1:p) {
		V[i] <- mean(xiljn[i,]*specs)
	}
	y <- do.call(g,as.list(V))
	len <- length(y)
	res <- matrix(nrow=B,ncol=len)
	cnames <- names(y)
		
	# Fit the linear autoregressive sieve process
	fit <- arima(x, order=c(ln,0,0))
	coeffs <- coef(fit)
	psi <- coeffs[1:ln]
	eps_hat <- residuals(fit)
	# Remove the first ln values and center on the mean
	eps_hat <- eps_hat[-(1:ln)]
	eps_tilde <- eps_hat - sum(eps_hat)/(n-ln)

	# Bootstrap
	for (j in 1:B) {
		# Draw epsilon values randomly with replacement
		eps_star <- sample(eps_tilde,n-ln,replace=TRUE)
		
		# Reconstruct the bootstrapped series
		xstar <- numeric(n)
		xstar[1:ln] <- x[1:ln]
		for (i in (ln+1):n) {
			xstar[i] <- sum(psi*xstar[(i-1):(i-ln)]) + eps_star[i-ln]
		}
		
		# Compute the periodograms of the xstar series
		P <- spec.pgram(ts(xstar,frequency=1),plot=FALSE,taper=0,fast=FALSE,detrend=FALSE)
		I_star_tilde <- P$spec/(2*pi)

		# Compute \hat{f} for the Fourier frequencies
		sigma2 <- sum(eps_star^2)/(n-ln)
		fhat <- numeric(n0)
		J <- 1:ln
		for (i in 1:n0) {
			f <- freqs[i]
			eijf <- complex(modulus=1,argument=-f*J)
			fhat[i] <- Mod(1-sum(psi*eijf))^(-2)*sigma2/(2*pi)
		}
		# Complete fhat by symmetry
		fhat0 <- Mod(1-sum(psi))^(-2)*sigma2/(2*pi)
		fhat <- c(rev(fhat),fhat0,fhat)

		# Compute \hat{q} for the Fourier frequencies
		QIf <- I_n/fhat
		qhat <- numeric(n0)
		for (i in 1:n0) {
			K <- smc[(n0-i+1):(3*n0-i+1)]
			qhat[i] <- 2*pi*mean(K*QIf)/h
		}
		I_star <- I_star_tilde*qhat
		V <- numeric(p)
		for (i in 1:p) {
			V[i] <- mean(xiljn[i,]*I_star)
		}
		res[j,] <- do.call(g,as.list(V))
	}
	
	if (len == 1) {
		res <- as.vector(res)
	} else if (!is.null(cnames)) {
		colnames(res) <- cnames
	}
	obj <- list(s=res,Tn=y)
	class(obj) <- "boodd"
	attr(obj,"kind") <- "aided"
	attr(obj,"func") <- g
    return(obj)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "bootglm(model,data,func,B,...)" --
 # 
 # Bootstrap of generalized linear model parameters.
 # 
 # ------------------------------------------------------------------------
 ##
bootglm <- function(model,data,func,B,...) {
	if (class(model)[1] != "glm" & class(model)[1] != "lm") {
		stop("'model' must be an object of class 'glm' or 'lm'")
	}
	
	# Test the value returned by func
	y <- func(data,...)
	if (!is.vector(y)) {
		stop("Function 'func' must return a vector")
	}
	len <- length(y)
	cnames <- names(y)
	res <- matrix(nrow=B,ncol=len)
	# Get the name of the response variable
	if (class(model)[1] == "glm") {
		repname <- all.vars(model$formula)[[1]]
	} else {
		repname <- all.vars(model$terms)[[1]]
	}
	
	# Generate B simulations
	sims <- simulate(model,nsim=B)
	for (i in 1:B) {
		data[[repname]] <- sims[[i]]
		res[i,] <- func(data,...)
	}	
	
	obj <- list(s=res,Tn=y)
	class(obj) <- "boodd"
	attr(obj,"kind") <- "bootglm"
	attr(obj,"func") <- func
    return(obj)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "fieldboot(arr,func,B,blocklens,method,...)" --
 # 
 # Bootstrap of discrete multidimensional random fields.
 # 
 # ------------------------------------------------------------------------
 ##
fieldboot <- function(arr,func,B,blocklens,...,method=c("movingblock","nonoverlapping","circular")) {
	# Test the value returned by func
	y <- func(arr,...)
	if (!is.vector(y)) {
		stop("Function 'func' must return a vector")
	}
	len <- length(y)
	cnames <- names(y)
	res <- matrix(nrow=B,ncol=len)
	
	# Data lengths n_1, n_2, ..., n_d
	dlens <- dim(arr)
	ndims <- length(dlens)
	if (ndims < 2) {
		stop("expected at least 2-dimensional array")
	}
	if (length(blocklens) == 1) {
		blocklens = rep(blocklens,ndims)
	}
	if (length(blocklens) != ndims) {
		stop("wrong number of block lengths")
	}
	if (any(blocklens > dlens)) {
		stop("block lengths must be less than data lengths in all dimensions")
	}
	if (any(blocklens <= 0)) {
		stop("block lengths must be positive")
	}
	method <- match.arg(method)
	# Block lengths b_1, b_2, ..., b_d
	blens <- blocklens
	last <- dlens-blens+1
	fulldims <- ((dlens-1)%/%blens +1)*blens
	
	# Array indices as string "1:n_1,1:n_2,...,1:n_d"
	indstr <- blockString(rep(1,ndims),dlens)
	
	# Build the tile indices (as strings) defining the covering
	# of the reconstructed array
	tiles <- blockTiles(blens,dlens)
	tnum <- length(tiles)

	# Initialize the new array (a slightly bigger array containing an 
	# entire number of tiles that will be trimmed later)
	narr <- array(0,dim=fulldims)
	origs <- integer(ndims)
	if (method[1] == "circular") {
		# Ensure circularity
		arr <- completeArray(arr,blens)	
	}
	
	for (i in 1:B) {
		for (j in 1:tnum) {
			tile <- tiles[j]
			
			# Draw a random block in the original array
			for (k in 1:ndims) {
				if (method[1] == "movingblock") {
					last <- dlens[k]-blens[k]+1
					origs[k] <- sample(1:last,1)
				} else if (method[1] == "nonoverlapping") {
					origs[k] <- sample(1+blens[k]*(0:(dlens[k]%/%blens[k]-1)),1)
				} else if (method[1] == "circular") {
					last <- dlens[k]
					origs[k] <- sample(1:last,1)
				}  
			}
			rndblock <- blockString(origs,blens)
			
			# Build the command string
			cmd <- paste0("narr[",tile,"] <- ","arr[",rndblock,"]")
			# Eval the command
			eval(parse(text=cmd))
		}
		
		# Apply the statistics to the rebuilt array (trimmed to the
		# original dimensions)
		resCmd <- paste0("res[",i,",] <- func(narr[",indstr,"],...)")
		eval(parse(text=resCmd))
	}
	
	if (len == 1) {
		res <- as.vector(res)
	} else if (!is.null(cnames)) {
		colnames(res) <- cnames
	}
	obj <- list(s=res,Tn=y)
	class(obj) <- "boodd"
	attr(obj,"kind") <- c("randomfield",method[1])
	attr(obj,"func") <- func
    return(obj)
}

