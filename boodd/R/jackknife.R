# ===========================================================================
# File: "jackknife.R"
#                        Created: 2020-10-20 12:18:43
#              Last modification: 2020-11-19 11:23:42
# Author: Bernard Desgraupes
# e-mail: <bernard.desgraupes@parisnanterre.fr>
# This file is part of the boodd project.
# ===========================================================================


## 
 # ------------------------------------------------------------------------
 # 
 # "jackVar(x,func)" --
 # 
 # Jackknife estimator for the variance of function 'func'.
 # 
 # ------------------------------------------------------------------------
 ##
jackVar <- function(x,func,...) {
	# Test the value returned by func
	Tn <- func(x,...)
	if (!is.vector(Tn)) {
		stop("Function 'func' must return a vector")
	}
	p <- length(Tn)
	Tnames <- names(Tn)
	if (is.matrix(x)) {
		n <- nrow(x)
	} else {
		n <- length(x)
	}
	if (n <= p) {stop("length of x must be greater than length of func(x)")}
	J <- matrix(nrow=n,ncol=p)
	for (i in 1:n) {
		if (is.matrix(x)) {
			Tni <- func(x[-i,],...)
		} else {
			Tni <- func(x[-i],...)
		}
		J[i,] <- n*Tn - (n-1)*Tni
	}
	muJ <- colMeans(J)
	V <- matrix(0,nrow=p,ncol=p)
	for (i in 1:n) {
		L <- J[i,]-muJ
		V <- V + L%*%t(L)
	}
	V <- V/(n*(n-p))
	if (!is.null(Tnames)) {
		rownames(V) <- Tnames
		colnames(V) <- Tnames
	}
	if (p == 1) {V <- as.vector(V)}
	return(V)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "jackVarBlock(x,func,blen,...)" --
 # 
 # Jackknife estimator by blocks for the variance of function 'func'.
 # 
 # ------------------------------------------------------------------------
 ##
jackVarBlock <- function(x,func,blen,...) {
	# Test the value returned by func
	Tn <- func(x,...)
	if (!is.vector(Tn)) {
		stop("Function 'func' must return a vector")
	}
	n <- length(x)
	if (n <= blen) {stop("blen must be less than length of x")}	
	bnum <- n%/%blen
	p <- length(Tn)
	Tnames <- names(Tn)	
	J <- matrix(nrow=bnum,ncol=p)
	for (i in 1:bnum) {
		# Indices of the block to remove
		bidx <- ((i-1)*blen+1):(i*blen)
		Tni <- func(x[-bidx],...)
		J[i,] <- (n*Tn - (n-blen)*Tni)/blen
	}
	muJ <- colMeans(J)
	V <- matrix(0,nrow=p,ncol=p)
	for (i in 1:bnum) {
		L <- J[i,]-muJ
		V <- V + L%*%t(L)
	}
	V <- V/(n*bnum)
	if (!is.null(Tnames)) {
		rownames(V) <- Tnames
		colnames(V) <- Tnames
	}
	if (p == 1) {V <- as.vector(V)}
	return(V)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "jackVarRegen(x,func,...,atom,small,s)" --
 # 
 # Jackknife estimator by regenerative blocks for the variance of function 'func'.
 # 
 # ------------------------------------------------------------------------
 ##
jackVarRegen <- function(x,func,...,atom,small=NULL,s=median(x)) {
	if (!missing(atom)) {
		# This is the atomic case
		res <- jackVarRegen.atom(x,func,atom=atom,...)
	} else {
		if (is.null(small)) {
			small <- findBestEpsilon(x,s)
		}
		res <- jackVarRegen.smallEnsemble(x,func,small,...)
	}
	return(res)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "jackVarRegen.atom(x,func,atom,...)" --
 # 
 # Regenerative jackknife estimator of the variance of function
 # 'func' for finite states Markov chains.
 # 
 # ------------------------------------------------------------------------
 ##
jackVarRegen.atom <- function(x,func,atom,...) {
	if (length(atom) != 1) {
		stop("atom must be a single value")
	}
	if ( !(atom %in% x) ) {
		stop("atom must be an element of x")
	}
	Tn <- func(x,...)
	if (!is.vector(Tn)) {
		stop("Function 'func' must return a vector")
	}
    n <- length(x)
	nTn <- n*Tn
	p <- length(Tn)
	Tnames <- names(Tn)	

	idx <- which(x==atom)
	li <- length(idx)
	starts <- idx+1
	lens <- diff(c(idx,n))
	# Drop the last block, either because it is incomplete or because it is
	# past the last position
	nb <- li - 1
	
	J <- matrix(nrow=nb,ncol=p)
	for (i in 1:nb) {
		start <- starts[i]
		len <- lens[i]
		# Indices of the block to remove
		bidx <- start:(start+len-1)
		Tni <- func(x[-bidx],...)
		J[i,] <- (nTn - (n-len)*Tni)/len
	}
	muJ <- colMeans(J)
	V <- matrix(0,nrow=p,ncol=p)
	for (i in 1:nb) {
		L <- J[i,]-muJ
		V <- V + L%*%t(L)
	}
	V <- V/(n*nb)
	if (!is.null(Tnames)) {
		rownames(V) <- Tnames
		colnames(V) <- Tnames
	}
	if (p == 1) {V <- as.vector(V)}		

	return(V)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "jackVarRegen.smallEnsemble(x,func,small,...)" --
 # 
 # Regenerative jackknife estimator of the variance of function
 # 'func' for homogeneous Markov chains.
 # 
 # ------------------------------------------------------------------------
 ##
jackVarRegen.smallEnsemble <- function(x,func,small,...) {
	if (class(small)[1] != "smallEnsemble") {
		stop("expecting an object of class 'smallEnsemble'.")
	}
	
	Tn <- func(x,...)
	if (!is.vector(Tn)) {
		stop("Function 'func' must return a vector")
	}
	n <- length(x)
	nTn <- n*Tn
	p <- length(Tn)
	Tnames <- names(Tn)	

	s <- small$s
	eps <- small$epsilon
	delta <- small$delta
	p_XiXip1 <- small$trans	

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
	
	J <- matrix(nrow=nb,ncol=p)
	for (i in 1:nb) {
		start <- S[i]
		len <- L[i]
		# Indices of the block to remove
		bidx <- start:(start+len-1)
		Tni <- func(x[-bidx],...)
		J[i,] <- (nTn - (n-len)*Tni)/len
	}
	muJ <- colMeans(J)
	V <- matrix(0,nrow=p,ncol=p)
	for (i in 1:nb) {
		L <- J[i,]-muJ
		V <- V + L%*%t(L)
	}
	V <- V/(n*nb)
	if (!is.null(Tnames)) {
		rownames(V) <- Tnames
		colnames(V) <- Tnames
	}
	if (p == 1) {V <- as.vector(V)}

	return(V)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "jackFunc(func)" --
 # 
 # Create a vector-valued function that calculates both the statistics
 # defined by 'func' and the estimated jackknife variance. Return a
 # function object.
 # 
 # Example of use:
 #     func <- function(x) {sum(abs(x))}
 #     jf <- jackFunc(func)
 #     boo <- boots(x,jf,10)
 # 
 # ------------------------------------------------------------------------
 ##
jackFunc <- function(func,...) {
	jf <- function(X) {
		jv <- jackVar(X,func,...)
		# diag() expects a matrix
		if (!is.matrix(jv)) {jv <- as.matrix(jv)}
		return(c(func(X,...),diag(jv)))
	}
	return(jf)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "jackFuncBlock(func,blen)" --
 # 
 # Create a vector-valued function that calculates both the statistics
 # defined by 'func' and the estimated jackknife variance by blocks. Return
 # a function object.
 # 
 # Example of use:
 #     func <- function(x) {sum(abs(x))}
 #     jfb <- jackFuncBlock(func,blen=10)
 #     boo <- boots(x,jfb,99)
 # 
 # ------------------------------------------------------------------------
 ##
jackFuncBlock <- function(func,blen=NULL,...) {
	jfb <- function(X) {
		if (is.null(blen)) {
			blen <- floor(length(X)^(1/3))
		}
		jv <- jackVarBlock(X,func,blen=blen,...)
		# diag() expects a matrix
		if (!is.matrix(jv)) {jv <- as.matrix(jv)}
		return(c(func(X,...),diag(jv)))
	}
	return(jfb)
}


## 
 # ------------------------------------------------------------------------
 # 
 # "jackVarField(arr,func,blocklens,...)" --
 # 
 # Jackknife estimator of the variance of function 'func' for random fields.
 # 
 # ------------------------------------------------------------------------
 ##
jackVarField <- function(arr,func,blocklens,...) {
	if (any(is.nan(arr))) {
		stop("array contains NaNs")
	}
	# Test the value returned by func
	Tn <- func(arr,...)
	if (!is.vector(Tn)) {
		stop("Function 'func' must return a vector")
	}
	p <- length(Tn)
	Tnames <- names(Tn)	
	dlens <- dim(arr)
	n <- prod(dlens)
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
	blens <- blocklens
	bsize <- prod(blens)

	# Build the tile indices as strings defining the covering
	# of the array (exceeding values are trimmed)
	tiles <- blockTiles(blens,dlens,canExceed=FALSE)
	tnum <- length(tiles)
	V <- as.vector(arr)

	J <- matrix(nrow=tnum,ncol=p)
	for (i in 1:tnum) {
		narr <- arr
		# Indices of the current tile as string
		tile <- tiles[i]
		# Replace values to remove by NaN
		cmd <- paste0("narr[",tile,"] <- NaN")
		eval(parse(text=cmd))
		# Indices of the items to remove		
		bidx <- which(is.nan(narr))
		Tni <- func(V[-bidx],...)
		J[i,] <- (n*Tn - (n-bsize)*Tni)/bsize
	}
	muJ <- colMeans(J)
	V <- matrix(0,nrow=p,ncol=p)
	for (i in 1:tnum) {
		L <- J[i,]-muJ
		V <- V + L%*%t(L)
	}
	V <- V/(n*tnum)
	if (!is.null(Tnames)) {
		rownames(V) <- Tnames
		colnames(V) <- Tnames
	}
	if (p == 1) {V <- as.vector(V)}
	return(V)
}
