###############################
##
## This file contains:
##
##	rlikfit()
##	print.rlikfit()
##	summary.rlikfit()
##	print.summary.rlikfit()
##  prof.rlik()

##' @name rlikfit
##' @title Likelihood Fit for Regoinalized Universal Kriging Model
##
##' @description Computes maximum likelihood estimates for parameters
##'			in a universal kriging model with regionally-varying
##'			covariance parameters.
##
##' @details The MLEs are obtained by maximizing the profile
##'		likelihood (see \code{\link{prof.rlik}}).
##'  By default, the opimization is done using the 'L-BFGS-B' method of
##'  \code{\link{optim}}. 
##'
##'  If provided, the initial (log) parameters should be given 
##'	 as a vector of length 3 times the number of distinct regions.
##'	 Values should follow the order of nugget,
##'	 sill, and range for each region in succession. 
##'  For example, if there are three regions the vector 
##'  c(tau1, sigma1, rho1, tau2, sigma2, rho2, tau3, sigma3, rho3) should 
##'	 be given. If unspecifed, default initial values of 0 (on the log scale) 
##'  are used.
##'  The default upper bound for the nugget and sill is log(10), and
##'  the default lower bound is log(0.0001). Default bounds for the range
##'  are log(0.001) and log(5000). These can be modified by providing 
##'  'upper' and 'lower' values in \code{optim.args}.
##'
##'  Regionally-varying coefficients for regression parameters
##'	 are not explicitly implemented. However, they can be obtained
##'	 by supplying the covariate matri \code{X} with suitable block
##'	 structure.  
##'
##'	 Setting \code{cov.model='iid'} for all regions is equivalent
##'	 to ordinary least squares regression, but this function is
##'	 is not currently optimized for this use case and will likely
##'	 be much slower and less precise than \code{\link{lm}}.
##
## Input:
##
##' @param y Dependent variable to be modeled
##' @param  X  A \eqn{n x d} matrix of kriging covariates
##'      for the mean component.
##' @param coords \eqn{N x 2} matrix of coordinates
##' @param reg.ind \eqn{N}-vector of distinct kriging regions.  
##'		Defaults to  \code{reg.ind = rep(1,N)}, 
##'		which yields traditional universal kriging.
##' @param cov.model Either a single character string or 
##'		a named vector providing the covariance model to be used.
##'     Valid values are "exp" and "iid". If only one value is provided,
##'     it is used for all regions.
##' @param  init.pars  Vector of initial starting values 
##'    for log-covariance parameters.
##' @param optim.args List of named arguments passed to \code{\link{optim}}. 
##
##' @author Joshua Keller, Paul Sampson
##' @export
##
##  Output:
##' @return An object of class 'rlikfit' containing:
##' \item{log.cov.pars}{maximized log-covariance parameters}
##'  \item{beta}{kriging regression covariate estimates}
##'  \item{max.log.lik}{maximized log-likelihood}
##'  \item{hess.pd}{ indicator of positive definite hessian, if hess=TRUE}
##'  \item{r}{ Number of regions}
#
##' @importFrom stats optim
rlikfit <- function(y, X, coords, reg.ind, cov.model="exp", init.pars=NULL, optim.args=NULL)
{
	if (missing(reg.ind)){
		reg.ind <- rep(1, length(y))
	}
	u.reg <- unique(reg.ind)
	r <- length(u.reg)
	ind.order <- order(reg.ind)

	if (length(cov.model)==1){
		cov.model <- rep(cov.model, times=r)
		names(cov.model) <- u.reg
	} 
	if (length(cov.model)!=r){
		stop("Provided 'cov.model' not equal to the number of regions.")
	}
	if (!all(names(cov.model) %in% u.reg)){
		stop("Names of 'cov.model' should match unique elements of 'reg.ind'.")
	}
	cov.model <- cov.model[unique(reg.ind[ind.order])]
	
	npars <- sum(c("exp"=3, "iid"=1)[match(cov.model, c("exp", "iid"))])
	if (is.null(init.pars)){
		init.pars <- rep(log(1), times= npars)
	}
	if (length(init.pars)!=npars){
		stop("init.pars has incorrect length.")
	}

	X <- X[ind.order,, drop=FALSE]
  	coords <- coords[ind.order,]
  	y <- y[ind.order]
  	reg.ind <- reg.ind[ind.order]
  	n.vec <- as.vector(by(reg.ind,reg.ind,length))

  	# Set default optim values.
  	optim.args.default <- list(method="L-BFGS-B", lower=rep(c(log(0.0001),log(0.0001),log(0.001)),length(n.vec)),upper=rep(c(log(10),log(10),log(5000)),length(n.vec)),hessian=TRUE)
	defaultargs <- which(!names(optim.args.default) %in% names(optim.args))
	optim.args <- c(optim.args, optim.args.default[defaultargs])

  opt <- do.call(stats::optim, args=c(list(par= init.pars, fn=prof.rlik,coords=coords, y=y, X=X, reg.ind=reg.ind, cov.model= cov.model), optim.args))
  	beta <- inf.beta(pars = opt$par, X = X, coords = coords, 
                   reg.ind = reg.ind, y  = y, cov.model=cov.model)
  	eigs <- eigen(opt$hessian, symmetric=T, only.values=T)$values
  	# Note that this is hessian of *negative* log likelihood, so
  	# positive def = negative definite hessian = maximum.
  	hess.pd <- sum(eigs>(length(eigs)*.Machine$double.eps))==length(eigs)
  	if (!hess.pd) warning("Hessian not negative definite.")
  out <- list('log.cov.pars'=opt$par,'beta'=beta,'max.log.lik'=-opt$value,'hess.pd'= hess.pd, 'r'=r, 'cov.model'=cov.model, 'opt'=opt)
  class(out) <- "rlikfit"
  out
}



##' @title Print details for class \code{rlikfit}
##' @description \code{\link[base:print]{print}} method for class \code{rlikfit}.
##' @param x object of class \code{rlikfit}
##' @param ... Ignored additional arguments.
##' @export
##' @family 'rliklfit methods'
print.rlikfit <- function(x, ...) {
	if(class(x)!="rlikfit"){
		stop("x must be of class 'rlikfit'.")
	}	
	cat("An object of class rlikfit.\n")
	cat("Estimated log covariance parameters are:\n")
	cov.pars <- matrix(NA, nrow=x$r, ncol=3)
	ii <- 0
	for (i in 1:x$r){
		cov.inds <- if(x$cov.model[i]=="exp") {ii + 1:3} else {ii + 1}
		cov.pars[i,1:length(cov.inds)] <- x$log.cov.pars[cov.inds]
		ii <- max(cov.inds)
	}
	pm <- matrix(cov.pars, ncol=3, byrow=TRUE, dimnames=list(paste0("Region ", 1:x$r), c("Tau", "Sigma", "Phi")))
	rownames(cov.pars) <- names(x$cov.model)
	colnames(cov.pars) <- c("Tau", "Sigma", "Phi")
	print(cov.pars)
} ##print.rlikfit()

##' @title Compute summary details for class \code{rlikfit}
##' @description \code{\link[base:print]{summary}} method for class \code{rlikfit}.
##' @param object object of class \code{rlikfit}
##' @param ... Ignored additional arguments.
##' @export
summary.rlikfit <- function(object, ...) {
	
	class(object) <- "summary.rlikfit"
	return(object)
}

##' @export
print.summary.rlikfit <- function(x, ...){
	cat("An object of class rlikfit.\n")
	cat("Estimated log covariance parameters are:\n")
	pm <- matrix(NA, nrow=x$r*3, ncol=2)	
	cov.ind.count <- sapply(ifelse(x$cov.model=="exp", 3, 1), function(w) 1:w)	
	cov.ind <- unlist(mapply("+", cov.ind.count, 3*(0:(x$r-1))))
	pm[cov.ind, 1] <- x$log.cov.pars
	rownames(pm) <- paste(rep(names(x$cov.model), each=3), rep(c("Tau", "Sigma", "Phi"), times=x$r))
	pm[cov.ind, 2] <- suppressWarnings(sqrt(diag(x$opt$hessian)))
	print(pm)
	if (any(is.nan(pm[,2]))) {
		cat("Warning: Hessian not positive definite, MLE\n")
		cat("         may not have been reached.\n")}
	cat("\nEstimated regression coefficients:\n")
	print(t(x$beta))
}





#=================================================================
# 
#
# likfit.obj: Object output from my.likfit()
# X.mon: Covariates (i.e., PLS scores) for monitors
# coords.mon: Monitor coordinates
# reg.mon: Vector of regions for monitors
# y: Monitor pollutant concentrations
# X.pred: Covariates (i.e., PLS scores) for prediction locations
# coords.pred: Prediction coordinates
# reg.pred: Vector of regions for prediction sites
#
# Returns: N x 2 matrix, where N is the number of prediction locations;
# 1st column contains predictions, 2nd column contains prediction variances
#====================================================================
#====================================================================
##' @name predict.rlikfit
##' @title Predict from a regionalized kriging model
##
##' @description Makes conditional expectations (kriging predictions) from
##' 	a regionalized kriging model and computes prediction variances.
##
## @details Details to be added.
##' @author Joshua Keller, Paul Sampson
##'
##' @param object An object of class \code{\link{rlikfit}} from which to compute predictions.
##' @param X.mon Kriging covariates for monitors
##' @param coords.mon Monitor coordinates
##' @param reg.mon Vector of regions for monitors
##' @param y Monitor pollutant concentrations
##' @param X.pred Covariates (i.e., PLS scores) for prediction locations
##' @param coords.pred Prediction coordinates
##' @param reg.pred Vector of regions for prediction sites
##' @param ... ignored additional arguments
##
##'	@export
##' @importFrom stats model.matrix
##'	@seealso rlikfit
predict.rlikfit <- function(object, X.mon, coords.mon, reg.mon, y, X.pred, coords.pred, reg.pred, ...)
{
  log.cov.pars <- object $log.cov.pars
  beta <- object $beta
  n1 <- nrow(coords.mon)
  n2 <- nrow(coords.pred)
  p.ind <- c(rep(1,n1),rep(2,n2))
#  X.all <- stats::model.matrix(~ -1+ rbind(X.mon,X.pred))
  X.all <- rbind(X.mon,X.pred)
  dim.X <- ncol(X.all)
  coords.all <- rbind(coords.mon,coords.pred)
  reg.all <- c(reg.mon,reg.pred)
  y.all <- c(y,rep(NA,n2))
  v.all <- rep(NA,n1+n2)
  pred.id <- 1:n2
  mp.ind <- c(rep(0,n1),pred.id)
  
  dat.all <- cbind(coords.all,X.all,reg.all,mp.ind,p.ind,y.all,v.all)
  c.o <- dat.all[order(reg.all),]
  uq.regs <- unique(dat.all$reg.all)
  n.regs <- length(uq.regs)
  
  n.pars <- ifelse(object$cov.model=="exp", 3, 1)
  cum.pars <- cumsum(n.pars)
  
  regional.new <- NULL
  for (j in 1:length(uq.regs))
  {
      regional.old <- regional.new
      regional.dat <- c.o[which(c.o$reg.all==uq.regs[j]),]
     if (j==1){
       regional.new <- rbind(regional.old,
       k.pred(regional.dat,dim.X,beta,log.cov.pars[1:(cum.pars[j])], cov.model=object$cov.model[j]))
     } else {
       regional.new <- rbind(regional.old,
       k.pred(regional.dat,dim.X,beta,log.cov.pars[(cum.pars[j-1]+1):(cum.pars[j])], cov.model=object$cov.model[j]))
     }

  }
  
  predictions <- regional.new[match(pred.id,regional.new$mp.ind),]
  out <- cbind('predictions'=predictions$y.all,'pred.var'=predictions$v.all)
  out
}



#====================
#PROFILE LIKELIHOOD
#====================
# Modified version by J. Keller
# That does explicit matrix multiplication and direct calculation of the
# log density
# Found to be ~50% faster than original prof.lik
##' @name prof.rlik
##' @title Profile Likelihood for Kriging
##
##' @description Evaluates profile likelihood function for (regionalized) universal kriging model.
##
##' @details 
##'		This function is primarily intended to be called
##'		via \code{\link{rlikfit}}, and thus has limited 
##'		checking of arguments.
##'
##' @author Paul Sampson, Josh Keller
##'
##' @param pars Log covariance parameters
##' @param X Kriging covariate model matrix
##'	@param coords Coordinates
##' @param reg.ind Region indicator
##' @param y Observed outcome values.
##' @param cov.model Covariance model for each region.
##' @param useSTpackage Logical. If the SpatioTemporal package is available, 
##'		uses block cholesky functions within that package for
##'		improved performance.
##' @return Scalar value giving the *negative* log profile likelihood
##'		for the kriging model.
##' @export
##' @importFrom stats mahalanobis
##' @seealso \link{rlikfit}
prof.rlik <- function(pars, X, coords, reg.ind, y, cov.model="exp", useSTpackage=TRUE)
{
 if (useSTpackage  && !requireNamespace("SpatioTemporal", quietly = TRUE)) {
		useSTpackage <- FALSE
		warning("Argument 'useSTpackage' is TRUE, but SpatioTemporal package not available. Setting to FALSE.")
  }
  Sig <- block.Sig(pars,coords,reg.ind, cov.model)
  beta <- inf.beta(pars=pars, X= X, coords, reg.ind, y, cov.model=cov.model, Sig=Sig, useSTpackage= useSTpackage)

  if(useSTpackage){
  	n.vec <- as.vector(by(reg.ind,reg.ind,length))
	Sig <- SpatioTemporal::makeCholBlock(Sig, n.blocks=length(n.vec), block.sizes=n.vec)
	Sig.inv <- SpatioTemporal:: invCholBlock(Sig, n.blocks=length(n.vec), block.sizes=n.vec)
  } else {
	Sig <- chol(Sig)
	Sig.inv <- chol2inv(Sig)
  }

 distval  <- stats::mahalanobis(x=as.numeric(y), center=X%*%beta, cov=Sig.inv, inverted=T)
 logdet <- 2*sum(log(diag(Sig)))
 lik <- -(length(y)*log(2 * pi) + logdet + distval)/2
 -lik
}

# prof.rlik.orig <- function(x, R, coords, reg.ind, y)
# {
	# l.pars <- x
	# m.mat <- R
  # Sig <- block.Sig(l.pars,coords,reg.ind)
  # beta <- inf.beta(l.pars, X= m.mat, coords, reg.ind, y, Sig=Sig, useSTpackage=TRUE)
 # lik <- dmvnorm(y,m.mat%*%beta,Sig,log=T)
  # -lik
# }


