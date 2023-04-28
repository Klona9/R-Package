# ===========================================================================
# File: "utils.R"
#                        Created: 2020-10-20 12:18:43
#              Last modification: 2020-11-19 10:25:08
# Author: Bernard Desgraupes
# e-mail: <bernard.desgraupes@parisnanterre.fr>
# This file is part of the boodd project.
# ===========================================================================


## 
 # ------------------------------------------------------------------------
 # 
 # "valueForPos(x,pos)" --
 # 
 # If pos is in [n,n+1[, calculate the linear interpolation between x[n] and
 # x[n+1]. The interpolated value is:
 #     x[n]+(x[n+1]-x[n])*(pos-n)
 # 
 # ------------------------------------------------------------------------
 ##
valueForPos <- function(x,pos) {
	len <- length(x)
	if (pos <= 1) {
		n <- 1
	} else if (pos >= len) {
		n <- len-1
	} else {
		n <- floor(pos)
	}
	l <- x[n]
	r <- x[n+1]
	res <- l+(r-l)*(pos-n)
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "findBestEpsilon(x,s)" --
 # 
 # Find the best radius epsilon for the truncatures in small ensembles. The
 # function returns a list with the following computed values:
 #     * optimal epsilon
 #     * corresponding lower bound delta
 #     * estimated transition densities of Markov chain p_n(X_i,X_{i+1})
 # 
 # ------------------------------------------------------------------------
 ##
findBestEpsilon <- function(x,s=median(x),plotIt=FALSE) {
	n <- length(x)-1
	y <- x-s
	q <- quantile(abs(y),0.9)
	eps <- seq(from=0,q,len=200)[-1]
	leps <- length(eps)
	
	delta <- numeric(leps)
	N_eps <- numeric(leps)
	h <- bandw1(x)

	# Estimator of the transition densities p_n(X_i,X_{i+1})
	p_XiXip1 <- fastNadaraya(x,h)
	if (min(p_XiXip1) == 0) {
		loc0 <- (p_XiXip1 == 0)
		# To avoid NaN when dividing by p_XiXip1
		p_XiXip1[loc0] <- 1
		p_XiXip1[loc0] <- min(p_XiXip1)/10
	} 
	
	# Kernel estimator of the transition. S is a subdivision of the
	# interval [s-q,s+q] (q is max(eps)).
	S <- seq(from=s-q,to=s+q,len=500)
	f <- naradamar(S,S,x,h)
	
	# Compute the estimation of the number of blocks for all the values of
	# epsilon
	for (j in 1:leps) {
		e <- eps[j]
		l_eps <- sum(S <= s-e)
		u_eps <- sum(S <= s+e)
		delta[j] <- min(f[l_eps:u_eps,l_eps:u_eps])
		# Indicator if pairs (X_i,X_{i+1}) located inside the interval [s-epsilon,s+epsilon]
		voisin <- (x[1:n]>=s-e)*(x[2:(n+1)]>=s-e)*(x[1:n]<=s+e)*(x[2:(n+1)]<=s+e)
		# Estimator of the number of regenerative blocks
		N_eps[j]=delta[j]*sum(voisin/p_XiXip1)
	}
	
	optidx <- which.max(N_eps)
	
	if (plotIt) {
		plot(eps,N_eps,type="l",xlab="epsilon",ylab=expression(N[eps]))
		abline(h=0,col="gray")
		abline(v=eps[optidx],col="red")
		text(eps[optidx],1,round(eps[optidx],4))
	}
	res <- list(s=s,epsilon=eps[optidx],delta=delta[optidx],trans=p_XiXip1)
	class(res) <- "smallEnsemble"
	return(res)
}
	

## 
 # ------------------------------------------------------------------------
 # 
 # "smallEnsemble(x,s=median(x),eps,delta)" --
 # 
 # Build manually an object of class 'smallEnsemble'. See also the function
 # findBestEpsilon which also computes the same kind of objects.
 # 
 # ------------------------------------------------------------------------
 ##
smallEnsemble <- function(s,eps,delta,trans) {
	res <- list(s=s,epsilon=eps,delta=delta,trans=trans)
	class(res) <- "smallEnsemble"
    return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "fastNadaraya(x,h)" --
 # 
 # Compute the Nadaraya-Watson estimator p_n(X_i,X_{i+1}) of the transition
 # density p(x,y).
 # 
 # ------------------------------------------------------------------------
 ##
fastNadaraya <- function(x,h) {
	n <- length(x)-1
	Xi <- x[1:n]
	Xip1 <- x[2:(n+1)]
	f <- numeric(n)
	for (i in 1:n) {
		f[i] <- mean(dnorm((x[i]-Xi)/h)*dnorm((x[i+1]-Xip1)/h))/mean(dnorm((x[i]-x)/h))/h
	}
	return(f)
}

# # See also in package stats
# # The Nadaraya–Watson kernel regression estimate.
# # ksmooth(x, y, kernel = c("box", "normal"), bandwidth = 0.5,
# #         range.x = range(x),
# #         n.points = max(100L, length(x)), x.points)


## 
 # ------------------------------------------------------------------------
 # 
 # "naradamar(Sx,Sy,x,h)" --
 # 
 # Compute the Naradamar kernel on a grid given by subdivisions Sx and Sy.
 # 
 # ------------------------------------------------------------------------
 ##
naradamar <- function(Sx,Sy,x,h) {
	n <- length(x)-1
	Xi <- x[1:n]
	Xip1 <- x[2:(n+1)]
	lx <- length(Sx)
	ly <- length(Sy)
	f <- matrix(0,lx,ly)
	for (k in 1:lx) {
		for (l in 1:ly) {
			f[k,l] <- mean(dnorm((Sx[k]-Xi)/h)*dnorm((Sy[l]-Xip1)/h))/mean(dnorm((Sx[k]-Xi)/h))/h
		}
	}
	return(f)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "bandw1(y)" --
 # 
 # Compute the bandwith for the Nadaraya-Watson and the Naradamar estimators.
 # 
 # ------------------------------------------------------------------------
 ##
bandw1 <- function(y,C=(2*pi*sqrt(2))^(-0.2)) {
	# IQR with type 4 corresponds to the linear interpolation
	# of the empirical cdf
	iqr <- IQR(y,type=4)
	n <- length(y)
	s <- sd(y)
	return(C*min(s,(iqr/1.34)/n^0.2))
}


## 
 # ------------------------------------------------------------------------
 # 
 # "smoothingCoefficients <- function(n0,h,kernel)" --
 # 
 # Return the quantities K(\lambda_{jn}/h) for j= -n,...,+n. The returned
 # vector has length 2*(n%/%2)+1 and is symmetrical. The ((n%/%2)+1)-th value
 # corresponds to frequency 0.
 # 
 # ------------------------------------------------------------------------
 ##
smoothingCoefficients <- function(n,h,kernel="normal") {
	N <- 2*(n%/%2)
	# Adjust the kernel name
	kernel <- match.arg(kernel,c("normal","box","epanechnikov"))
	kernel <- substr(kernel,1,3)
	# Compute the coefficients
	if (kernel == "nor") {
		Q <- 2*pi*(0:N)/(h*n)
		K <- dnorm(Q)
	} else {
		K <- numeric(N+1)
		jm <- floor(h*n/(2*pi))
		idx <- 0:jm
		Q <- 2*pi*idx/(h*n)
		if (kernel == "epa") {
			# Indices are 1-based
			K[idx+1] <- 0.75*(1-Q^2)
		} else {
			K[idx+1] <- rep(0.5,length(Q))
		}
	}
	# Complete by symmetry
	K <- c(rev(K[-1]),K)
	return(K)
}


### 
 # ------------------------------------------------------------------------
 # 
 # "qVar(x,alpha,hn,kernel)" --
 # 
 # Quantile variance. If qa is the quantile for alpha, the quantile
 # variance is given by:
 #         alpha(1-alpha)
 #     v = --------------
 #         fhat_hn(qa)^2
 # where
 #     fhat_hn(x) = 1/(n hn) \sum_1^n K( (x-X_i)/hn )
 # 
 # The argument 'alpha' can be a scalar or a vector of probabilities.
 # 
 # ------------------------------------------------------------------------
 ##
qVar <- function(x,alpha,hn=NULL,kernel=c("gaussian","epanechnikov","rectangular")) {
	if (!is.vector(x)) {
		stop("x must be a vector")
	}
	n <- length(x)
	if (is.null(hn)) {
		# Get the cross-validation bandwidth estimation
		s <- sd(x)
		hn = bw.ucv(x,lower=s*0.000001,upper=s*2)
		# The estimate provided by unbiased cross-validation is of order
		# Cn^(-1/5). We multiply by n^(-2/15) because we need it to be of
		# order n^(-1/3).
		hn <- hn*n^(-2/15)
	}
	kernel <- match.arg(kernel)
	Qs <- quantile(x,alpha)	
	lq <- length(Qs)
	v <- numeric(lq)
	for (i in 1:lq) {
		cr <- (Qs[i]-x)/hn
		K <- switch(kernel,
			gaussian = dnorm(cr,sd=hn),
			epanechnikov = ifelse(abs(x)<1, .75*(1-cr^2), 0),
			rectangular = ifelse(abs(absx)<1, .5, 0)
		)
		fhat <- mean(K)/hn
		v[i] <- alpha[i]*(1-alpha[i])/fhat^2	
	}
	return(v)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "blockString(origs,blens)" --
 # 
 # Build a string representing a block extracted from an array. 'origs' is
 # the vector of start indices on all the dimensions. 'blens' is the vector
 # of block lengths.
 # 
 # Example:
 #     > blockString(c(4,9,5),c(2,3,5))
 #     [1] "4:5,9:11,5:9"
 # 
 # ------------------------------------------------------------------------
 ##
blockString <- function(origs,blens) {
	ndims <- length(blens)

	substr <- character(0)
	for (i in 1:ndims) {
		substr <- c(substr,paste0(as.character(origs[i]),":",as.character(origs[i]+blens[i]-1)))
	}
	str <- paste0(substr,collapse=",")
	return(str)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "blockTiles(blens,dlens,canExceed)" --
 # 
 # Build all the tile indices covering the array. 'blens' are the block
 # lengths and 'dlens' the data lengths. Return the indices as a vector of
 # strings. The 'canExceed' argument has the same meaning as with the
 # segmentStrings function.
 # 
 # Example:
 #     > blockTiles(c(3,4,5),c(9,10,10))
 #     "1:3,1:4,1:5"   "4:6,1:4,1:5"   "7:9,1:4,1:5"  
 #     "1:3,5:8,1:5"   "4:6,5:8,1:5"   "7:9,5:8,1:5"  
 #     "1:3,9:12,1:5"  "4:6,9:12,1:5"  "7:9,9:12,1:5" 
 #     "1:3,1:4,6:10"  "4:6,1:4,6:10"  "7:9,1:4,6:10" 
 #     "1:3,5:8,6:10"  "4:6,5:8,6:10"  "7:9,5:8,6:10" 
 #     "1:3,9:12,6:10" "4:6,9:12,6:10" "7:9,9:12,6:10"
 #     > blockTiles(c(3,4,5),c(9,10,10),FALSE)
 #     "1:3,1:4,1:5"  "4:6,1:4,1:5"  "7:9,1:4,1:5" 
 #     "1:3,5:8,1:5"  "4:6,5:8,1:5"  "7:9,5:8,1:5" 
 #     "1:3,1:4,6:10" "4:6,1:4,6:10" "7:9,1:4,6:10"
 #     "1:3,5:8,6:10" "4:6,5:8,6:10" "7:9,5:8,6:10"
 # 
 # ------------------------------------------------------------------------
 ##
blockTiles <- function(blens,dlens,canExceed=TRUE) {
	bl <- length(blens)
	dl <- length(dlens)
	if (bl < 2) {
		stop("expected at least 2 dimensions")
	}
	if (bl != dl) {
		stop("'blens' and 'dlens' must have same length")
	}
	out <- outer(segmentStrings(blens[1],dlens[1],canExceed),segmentStrings(blens[2],dlens[2],canExceed),FUN="paste",sep=",")
	if (bl > 2) {
		for (i in 3:bl) {
			out <- outer(out,segmentStrings(blens[i],dlens[i],canExceed),FUN="paste",sep=",")
		}
	}
	return(as.vector(out))
}


## 
 # ------------------------------------------------------------------------
 # 
 # "segmentStrings(bl,dl,canExceed)" --
 # 
 # Build a vector of strings corresponding to the successive intervals along a 
 # direction. 'bl' is the block length and 'dl' the data length. The last 
 # range may exceed the data length (it will be trimmed in the fieldboot 
 # function) unless the 'canExceed' argument is set to FALSE.
 # 
 # Example:
 #     > segmentStrings(3,10)
 #     [1] "1:3"   "4:6"   "7:9"   "10:12"
 #     > segmentStrings(3,10,canExceed=FALSE)
 #     [1] "1:3"   "4:6"   "7:9"
 # 
 # ------------------------------------------------------------------------
 ##
segmentStrings <- function(bl,dl,canExceed=TRUE) {
	nb <- (dl-1)%/%bl+1
	if (!canExceed & dl%%bl > 0) {
		nb <- nb-1
	}
	str <- character(0)
	for (i in 1:nb) {
		str <- c(str,paste0(as.character(bl*(i-1)+1),":",as.character(bl*i)))
	}
	return(str)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "completeArray <- function(arr,blens)" --
 # 
 # Extend an array to handle circular blocks in random fields bootstrap. 
 # The function returns an array in which all the dimensions have been 
 # completed in order to be a multiple of the corresponding block length.
 # 
 # For instance, if 'arr' has dim 2x3x4 and if the block lengths are 2,2,3 
 # then the returned array has dim 2x4x6.
 # 
 # ------------------------------------------------------------------------
 ##
completeArray <- function(arr,blens) {
	ndim <- length(blens)
	perm <- seq(along=blens)
	perm <- c(perm[2:ndim],perm[1])
	perm <- c(2:ndim,1)

	for (i in 1:ndim) {
		arr <- completeLastDim(arr,blens[ndim])
		blens <- c(blens[2:ndim],blens[1])
		arr <- aperm(arr,perm)
	}
	return(arr)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "completeLastDim(arr,bl)" --
 # 
 # Complete the last dimension of an array by circularity to make it a
 # multiple of 'bl' (block length). 
 # 
 # ------------------------------------------------------------------------
 ##
completeLastDim <- function(arr,bl) {
	dims <- dim(arr)
	ndim <- length(dims)
# 	dl <- dims[ndim]
# 	compl <- bl-dl%%bl
	compl <- bl-1
	len <- prod(dims[-ndim])
	v <- as.vector(arr)
	v <- c(v,v[1:(compl*len)])
	newdims <- c(dims[-ndim],dims[ndim]+compl)
	narr <- array(v,dim=newdims)
	return(narr)
}

