# ===========================================================================
# File: "periodic.R"
#                        Created: 2018-07-31 11:21:48
#              Last modification: 2018-10-25 16:08:32
# Author: Bernard Desgraupes
# e-mail: <bernard.desgraupes@parisnanterre.fr>
# This file is part of the boodd project.
# ===========================================================================


seasonalMean <- function(x,d,...) UseMethod("seasonalMean")
seasonalVar <- function(x,d,...) UseMethod("seasonalVar")
seasonalACF <- function(x,tau,d,...) UseMethod("seasonalACF")
meanCoeff <- function(x,freq,...) UseMethod("meanCoeff")
acfCoeff <- function(x,tau,freq,...) UseMethod("acfCoeff")


## 
 # ------------------------------------------------------------------------
 # 
 # "seasonalMean.default(x,d,...)" --
 # 
 # Estimators of the seasonal means of a periodically correlated time series.
 #     \widehat{\mu}_s = \frac{1}{w_s} \sum_{i=0}^{w_s-1} X_{s+id}
 # Note that the seasons are always numbered from 1, i-e we do not take
 # into account the start of the series.
 # 
 # ------------------------------------------------------------------------
 ##
seasonalMean.default <- function(x,d,...) {
	x <- as.vector(x)
	n <- length(x)
	if (n < d) {
		d <- n
	} 
	F <- factor(rep(1:d,length.out=n))
	res <- tapply(x,F,mean)
	return(as.vector(res))
}


## 
 # ------------------------------------------------------------------------
 # 
 # "seasonalMean.ts(x,d,...)" --
 # 
 # Estimators of the seasonal means for a time series object. Argument 'd'
 # is optional and defaults to the frequency of the series.
 # 
 # ------------------------------------------------------------------------
 ##
seasonalMean.ts <- function(x,d=frequency(x),...) {
	return(seasonalMean.default(x,d,...))
}


## 
 # ------------------------------------------------------------------------
 # 
 # "seasonalMean.embb(x,d,...)" --
 # 
 # Estimators of the seasonal means of a bootstrap sample obtained by EMBB.
 #     \widehat{\mu}_s^{*}=\frac{1}{v^*_s}\sum_{i\in TI_s} X_{i}^{*}
 # 
 # ------------------------------------------------------------------------
 ##
seasonalMean.embb <- function(x,d,...) {
	if (missing(d)) {stop("period 'd' is required")}
	res <- numeric(d)
	S <- ((x[,2]-1)%%d) + 1
	for (s in 1:d) {
		I <- which(S==s)
		if (length(I)>0) {
			res[s] <- mean(x[I])
		} else {
			res[s] <- 0
		}
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "seasonalVar.default(x,d,...)" --
 # 
 # Estimators of the seasonal variances of a periodically correlated time series.
 #     \widehat{\sigma}_s^2 = \frac{1}{w_s} \sum_{i=0}^{w_s-1}
 #                                    ( X_{s+id}-\widehat{\mu}_s )^2
 # 
 # ------------------------------------------------------------------------
 ##
seasonalVar.default <- function(x,d,...) {
	muhat <- seasonalMean(x,d)
	x <- as.vector(x)
	n <- length(x)
	if (n<d) {
		d <- n
	} 
	res <- numeric(d)
	for (s in 1:d) {
		# Number of periodic indices s, s+d, s+2*d, ... up to n
		# This is E((n-s)/d) 
		nd <- floor((n-s)/d) 
		res[s] <- mean( (x[s+(0:nd)*d]-muhat[s])^2 )
	}		
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "seasonalVar.ts(x,d,...)" --
 # 
 # Estimators of the seasonal variances for a time series object. Argument 'd'
 # is optional and defaults to the frequency of the series.
 # 
 # ------------------------------------------------------------------------
 ##
seasonalVar.ts <- function(x,d=frequency(x),...) {
	return(seasonalVar.default(x,d,...))
}


## 
 # ------------------------------------------------------------------------
 # 
 # "seasonalVar.embb(x,d,...)" --
 # 
 # Estimators of the seasonal variances of a bootstrap sample obtained by
 # EMBB.
 #     \widehat{\sigma}_s^{*2}=\frac{1}{v_s^*}\sum_{i\in TI_s} ( X^{*}_i-\widehat{\mu}_s)^2
 # 
 # ------------------------------------------------------------------------
 ##
seasonalVar.embb <- function(x,d,...) {
	if (missing(d)) {stop("period 'd' is required")}
	X <- as.vector(x[,1])
	S <- ((x[,2]-1)%%d) + 1
	n <- length(X)
	res <- numeric(d)
	muhat <- seasonalMean(x,d)
	for (s in 1:d) {
		I <- which(S==s)
		if (length(I)>0) {
			res[s] <- mean( (X[I]-muhat[s])^2 )
		} else {
			res[s] <- 0
		}
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "seasonalACF.default(x,tau,d,...)" --
 # 
 # Estimators of the seasonal autocovariances of a periodically correlated
 # time series.
 #     \widehat{B}(s,\tau) = \frac{1}{w} \sum_{k=0}^{w-1}
 #                        ( X_{s+kd}-\widehat{\mu}_s ) ( X_{s+\tau+kd}-\widehat{\mu}_{<s+\tau>} )
 # 
 # ------------------------------------------------------------------------
 ##
seasonalACF.default <- function(x,tau,d,...) {
	if (!is.vector(tau) || !(round(tau)==tau) || any(tau<0)) {
		stop("argument 'tau' must be a vector of positive integers")
	}
	muhat <- seasonalMean(x,d)
	x <- as.vector(x)
	n <- length(x)
	w <- n%/%d
	lh <- length(tau)	
	res <- matrix(nrow=lh,ncol=d)
	base <- (0:(w-1))*d
	
	for (i in 1:lh) {
		h <- tau[i]
		# Elements of the sum for which $s+\tau+k*d>n$ are replaced by 0
		mxi <- w*d+h
		if (mxi > n) {
			x <- c(x,rep(0,mxi-n))
		}
		for (s in 1:d) {
			mus <- muhat[s]
			mush <- muhat[(s+h-1)%%d+1]
			res[i,s] <- mean( (x[base+s]-mus) * (x[base+s+h]-mush) )
		}
	}
	if (lh == 1) {
		return(as.vector(res))
	} else {
		rownames(res) <- tau
		colnames(res) <- 1:d
		return(res)
	}	
}


## 
 # ------------------------------------------------------------------------
 # 
 # "seasonalACF.ts(x,tau,d,...)" --
 # 
 # Estimators of the seasonal autocovariances for a time series object.
 # Argument 'd' is optional and defaults to the frequency of the series.
 # 
 # ------------------------------------------------------------------------
 ##
seasonalACF.ts <- function(x,tau,d=frequency(x),...) {
	return(seasonalACF.default(x,tau,d,...))
}


## 
 # ------------------------------------------------------------------------
 # 
 # "seasonalACF.embb(x,tau,d,...)" --
 # 
 # Estimators of the seasonal autocovariances of a bootstrap sample
 # obtained by EMBB.
 #     \widehat{B}^{*}(s,\tau)=\frac{1}{v_s^*}\sum_{t\in TI_s} 
 # 			          (X^{*}_{t}-\widehat{\mu}_sd) 
 # 			          (X^{*}_{t+\tau}-\widehat{\mu}_{<s+\tau>}d)
 # 
 # ------------------------------------------------------------------------
 ##
seasonalACF.embb <- function(x,tau,d,...) {
	if (missing(d)) {stop("period 'd' is required")}
	if (!is.vector(tau) || !(round(tau)==tau) || any(tau<0)) {
		stop("argument 'tau' must be a vector of positive integers")
	}
	X <- as.vector(x[,1])
	S <- ((x[,2]-1)%%d) + 1
	n <- length(X)
	w <- n%/%d
	lh <- length(tau)	
	res <- matrix(nrow=lh,ncol=d)
	muhat <- seasonalMean(x,d)
	
	for (i in 1:lh) {
		h <- tau[i]
		# Elements of the sum for which $s+\tau+k*d>n$ are replaced by 0
		len <- length(X)
		if (len < n+h) {
			X <- c(X,rep(0,n+h-len))
		}
		for (s in 1:d) {
			I <- which(S==s)
			if (length(I)>0) {
				mus <- muhat[s]
				mush <- muhat[(s+h-1)%%d+1]
				res[i,s] <- mean( (X[I]-mus) * (X[I+h]-mush) )
			} else {
				res[s] <- 0
			}
		}
	}
	if (lh == 1) {
		return(as.vector(res))
	} else {
		rownames(res) <- tau
		colnames(res) <- 1:d
		return(res)
	}	
}


## 
 # ------------------------------------------------------------------------
 # 
 # "meanCoeff.default(x,d,freq,...)" --
 # 
 # Estimators of the Fourier coefficients of the mean of a periodically
 # correlated time series.
 #     \widehat{b}(\gamma) = \frac{1}{n}\sum_{t=1}^n X_t e^{-i\gamma t}
 # 
 # ------------------------------------------------------------------------
 ##
meanCoeff.default <- function(x,d,freq=NULL,...) {
	if (missing(d)) {stop("period 'd' is required")}
	if (is.null(freq)) {
		freq <- 2*pi*(0:(d-1))/d
	} 	
	x <- as.vector(x)
	lf <- length(freq)	
	n <- length(x)
	t <- 1:n
	res <- numeric(lf)
	
	for (i in 1:lf) {
		f <- freq[i]
		eit <- complex(modulus=1,argument=-f*t)
		res[i] <- mean(x*eit)
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "meanCoeff.ts(x,d,freq,...)" --
 # 
 # Estimators of the Fourier coefficients of the mean for a time series object.
 # Argument 'd' is optional and defaults to the frequency of the series.
 # 
 # ------------------------------------------------------------------------
 ##
meanCoeff.ts <- function(x,d=frequency(x),freq=NULL,...) {
	return(meanCoeff.default(x,d,freq,...))
}


## 
 # ------------------------------------------------------------------------
 # 
 # "meanCoeff.embb(x,freq,...)" --
 # 
 # Estimators of the Fourier coefficients of the mean of a bootstrap sample
 # obtained by EMBB.
 #     \widehat{b}^*_n(\gamma)=\frac{1}{n}\sum_{t=1}^{n}X^*_t \exp(-i\gamma t^*)
 # 
 # ------------------------------------------------------------------------
 ##
meanCoeff.embb <- function(x,freq,...) {
	X <- as.vector(x[,1])
	t <- as.vector(x[,2])
	lf <- length(freq)	
	res <- numeric(lf)	
	for (i in 1:lf) {
		f <- freq[i]
		eit <- complex(modulus=1,argument=-f*t)
		res[i] <- mean(X*eit)
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "acfCoeff.default(x,tau,d,freq,...)" --
 # 
 # Estimators of the Fourier coefficients of the autocovariance of a
 # periodically correlated time series.
 #     \widehat{a}(\lambda,\tau) = \frac{1}{n}\sum_{t=1-\min\{\tau,0\}}^{n-\max\{\tau,0\}}
 # 			 (X_{t+\tau}-\widehat{\mu}_n(t+\tau))
 # 			 (X_{t}-\widehat{\mu}_n(t))e^{-i\lambda t}
 # 
 # ------------------------------------------------------------------------
 ##
acfCoeff.default <- function(x,tau,d,freq=NULL,...) {
	if (!is.vector(tau) || !(round(tau)==tau)) {
		stop("argument 'tau' must be a vector of integers")
	}
	if (is.null(freq)) {
		freq <- 2*pi*(0:(d-1))/d
	}
	lf <- length(freq)	
	n <- length(x)
	lh <- length(tau)	
	res <- matrix(nrow=lh,ncol=lf)
	
	# Calculate the mean estimator
	muhat <- numeric(n)
	bhat <- meanCoeff(x,freq,d)
	for (j in 1:n) {
		muhat[j] <- sum(bhat*complex(modulus=1,argument=-freq*j))
	}
	
	for (i in 1:lh) {
		h <- tau[i]
		if (h >= 0) {
			t <- 1:(n-h)
		} else {
			t <- (1-h):n
		}
		# Calculate the Fourier coefficients estimator of the acf for lag h
		for (j in 1:lf) {
			f <- freq[j]
			eit <- complex(modulus=1,argument=-f*t)
			res[i,j] <- mean((x[t+h]-muhat[t+h])*(x[t]-muhat[t])*eit)
		}
	}

	if (lh == 1) {
		return(as.vector(res))
	} else {
		rownames(res) <- tau
		colnames(res) <- paste("f",1:lf,sep="")
		return(res)
	}	
}


## 
 # ------------------------------------------------------------------------
 # 
 # "acfCoeff.ts(x,tau,d,freq,...)" --
 # 
 # Estimators of the Fourier coefficients of the autocovariance for a time
 # series object. Argument 'd' is optional and defaults to the frequency of
 # the series.
 # 
 # ------------------------------------------------------------------------
 ##
acfCoeff.ts <- function(x,tau,d=frequency(x),freq=NULL,...) {
	return(acfCoeff.default(x,tau,d,freq,...))
}


## 
 # ------------------------------------------------------------------------
 # 
 # "acfCoeff.embb(x,tau,freq,...)" --
 # 
 # Estimators of the Fourier coefficients of the autocovariance of a
 # bootstrap sample obtained by EMBB.
 # For the sake of simplicity we assume that E[X_t]\equiv 0.
 #     \widehat{a}^*_n (\lambda,\tau)=\frac{1}{n}\sum_{t=1}^{n-\tau}X^*_t X^*_{t+\tau}\exp (-i\lambda t^*)
 # 
 # ------------------------------------------------------------------------
 ##
acfCoeff.embb <- function(x,tau,freq,...) {
	if (!is.vector(tau) || !(round(tau)==tau)) {
		stop("argument 'tau' must be a vector of integers")
	}
	X <- as.vector(x[,1])
	t <- as.vector(x[,2])
	n <- length(X)
	lf <- length(freq)	
	lh <- length(tau)	
	res <- matrix(nrow=lh,ncol=lf)
	
	for (i in 1:lh) {
		h <- tau[i]
		if (h >= 0) {
			t <- 1:(n-h)
		} else {
			t <- (1-h):n
		}
		# Calculate the Fourier coefficients estimator of the acf for lag h
		for (j in 1:lf) {
			f <- freq[j]
			eit <- complex(modulus=1,argument=-f*t)
			res[i,j] <- mean(x[t]*x[t+h]*eit)
		}
	}

	if (lh == 1) {
		return(as.vector(res))
	} else {
		rownames(res) <- tau
		colnames(res) <- paste("f",1:lf,sep="")
		return(res)
	}	
}


## 
 # ------------------------------------------------------------------------
 # 
 # "embb.sample(x,length.block,method)" --
 # 
 # Extended moving block bootstrap (EMBB) for periodically correlated time
 # series (PC) and for almost periodically correlated time series (APC).
 # EMBB exists in 2 versions: standard (like MBB) and circular (like
 # circular MBB). Method names maybe abbreviated.
 # Note that, even for circular approach, we do not need to take care of
 # having an integer number of periods. We just take the rest of blocks
 # from the beginning of the sample, not thinking of periodicity.
 # 
 # ------------------------------------------------------------------------
 ##
embb.sample <- function(x,length.block,method=c("movingblock","circular")) {
	if (!is.ts(x)) {
		x <- ts(x)
	}
	if (length.block >= length(x)) {
		stop("the block length must be less than the size of the time series")
	}
	method <- match.arg(method)

	if (method[1] == "movingblock") {
		res <- moving.embb(x,length.block)
	} else if (method[1] == "circular") {
		res <- circular.embb(x,length.block)
	} 
	class(res) <- c("embb","matrix")
    return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "moving.embb(x,length.block)" --
 # 
 # Standard extended moving block method. The function returns a bivariate
 # series as a matrix whose first column is the bootstrapped x and second
 # column contains the original time indices of the chosen observations.
 # 
 # ------------------------------------------------------------------------
 ##
moving.embb <- function(x,length.block) {	
	n <- length(x)
	b <- length.block
	I <- 1:n
	# Number of blocks
	q <- n-b+1
	# Number of sections
	s <- ((n-1) %/% b)+1
	# Bivariate bootstrapped series
	res <- matrix(nrow=s*b,ncol=2)
	
	# Bootstrap the blocks (represented by their starting index)
	bb <- sample(1:q,size=s,replace=TRUE)
	pos <- 1
	# Rebuild a new series with the chosen blocks
	for (idx in bb) {
		res[pos:(pos+b-1),1] <- x[idx:(idx+b-1)]
		res[pos:(pos+b-1),2] <- I[idx:(idx+b-1)]
		pos <- pos+b
	}
	return(res[I,])
}


## 
 # ------------------------------------------------------------------------
 # 
 # "circular.embb(x,length.block)" --
 # 
 # Standard extended moving block method. The function returns a bivariate
 # series as a matrix whose first column is the bootstrapped x and second
 # column contains the original time indices of the chosen observations.
 # 
 # ------------------------------------------------------------------------
 ##
circular.embb <- function(x,length.block) {	
	n <- length(x)
	b <- length.block
	I <- 1:n
	# Ensure circularity
	x <- c(x,x[1:(b-1)])
	I <- c(I,1:(b-1))

	# Number of sections
	s <- ((n-1) %/% b)+1
	# Bivariate bootstrapped series
	res <- matrix(nrow=s*b,ncol=2)
	
	# Bootstrap the blocks (represented by their starting index)
	bb <- sample(1:n,size=s,replace=TRUE)
	pos <- 1
	# Rebuild a new series with the chosen blocks
	for (idx in bb) {
		res[pos:(pos+b-1),1] <- x[idx:(idx+b-1)]
		res[pos:(pos+b-1),2] <- I[idx:(idx+b-1)]
		pos <- pos+b
	}
	return(res[1:n,])
}


