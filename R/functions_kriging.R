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
##' @description Maximizes a kriging likelihood that can include a regionalized variogram and regionalized coefficients.
##
##' @details The maximization is done via profile likelihood, using the function \code{\link{prof.rlik}}.
##'
##'  The default upper bound for the nugget and sill is log(10), and the default lower bound is log(0.0001). Default bounds for the range are log(0.001) and log(5000). FINISH
#
#  Input:
##' @param  pars  Vector of initial starting values for log-covariance parameters. 
##'  Should have length \eqn{3r}, where	
##'  where \eqn{r} is the number of 
##'  distinct regions. Values should follow the order of nugget (\eqn{\tau}),
##'	 sill (\eqn{\sigma}), and range (\eqn{\phi}). 
##'  e.g. c(tau1,sigma1,rho1,tau2,sigma2,rho2,tau3,sigma3,rho3)
##' @param y Dependent variable to be modeled
##' @param  X  \eqn{n x d} matrix of kriging covariates for the mean component
##' @param coords \eqn{N x 2} matrix of coordinates
##' @param reg.ind \eqn{N}-vector of distinct kriging regions.  Defaults to  \code{reg.ind = rep(1,N)}, 
##'		which yields traditional universal kriging.
##' @param hess Logical; whether or not to compute hessian and test whether it is positive-definite.
## 	A positive-definite hessian is a good indicator that true log-likelihood max has been found
##' @param optim.args List of arguments passed to \code{\link{optim}}. See details.
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
rlikfit <- function(pars, y, X, coords, reg.ind,  hess=TRUE, optim.args=NULL)
{
  X <- X[order(reg.ind),, drop=FALSE]
  coords <- coords[order(reg.ind),]
  y <- y[order(reg.ind)]
  reg.ind <- reg.ind[order(reg.ind)]
  n.vec <- as.vector(by(reg.ind,reg.ind,length))
  r <- length(n.vec)
  
  optim.args.default <- list(method="L-BFGS-B",lower=rep(c(log(0.0001),log(0.0001),log(0.001)),length(n.vec)),upper=rep(c(log(10),log(10),log(5000)),length(n.vec)),hessian=TRUE)
defaultargs <- which(!names(optim.args.default) %in% names(optim.args))
optim.args <- c(optim.args, optim.args.default[defaultargs])

  opt <- do.call(stats::optim, args=c(list(par=pars, fn=prof.rlik,coords=coords, y=y, R=X, reg.ind=reg.ind), optim.args))
  beta <- inf.beta(l.pars = opt$par, X = X, coords = coords, 
                   reg.ind = reg.ind, y  = y)
  # Check PD of Hessian
  if (hess) {
  	eigs <- eigen(opt$hessian, symmetric=T, only.values=T)$values
  	hess.pd <- sum(eigs>(length(eigs)*.Machine$double.eps))==length(eigs) 
  }
  if (hess)
   out <- list('log.cov.pars'=opt$par,'beta'=beta,'max.log.lik'=-opt$value,'hess.pd'= hess.pd, 'r'=r)
 if (!hess) out <- list('log.cov.pars'=opt$par,'beta'=beta,'max.log.lik'=-opt$value, 'r'=r)
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
	pm <- matrix(x$log.cov.pars, ncol=3, byrow=F, dimnames=list(paste0("Region ", 1:x$r), c("Tau", "Sigma", "Phi")))
	print(pm)
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
	pm <- matrix(x$log.cov.pars, ncol=1)
	rownames(pm) <-  paste0("Region ", outer(1: x $r, c("Tau", "Sigma", "Phi"), paste))
	print(pm)
}






#=======================================================================
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
#=======================================================================
#=======================================================================
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
  X.all <- stats::model.matrix(~ -1+ rbind(X.mon,X.pred))
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
  
  regional.new <- NULL
  for (j in 1:length(uq.regs))
  {
      regional.old <- regional.new
      regional.dat <- c.o[which(c.o$reg.all==uq.regs[j]),]
      regional.new <- rbind(regional.old,
                            k.pred(regional.dat,dim.X,beta,log.cov.pars[(3*j-2):(3*j)]))
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
##' @title Evaluates profile likelihood function for (regionalized) universal kriging.
##
##' @description To do.
##
##' @details If the SpatioTemporal package is available, uses block cholesky
##'			functions within that package.
##' @author Paul Sampson, Josh Keller
##'
##' @param x Log covariance parameters
##' @param R Kriging covariate model matrix
##'	@param coords Coordinates
##' @param reg.ind Region indicator
##' @param y Observed outcome values.
##' @param useSTpackage See details.
##
##' @export
##' @importFrom stats mahalanobis
##' @seealso \link{rlikfit}
prof.rlik <- function(x, R, coords, reg.ind, y, useSTpackage=TRUE)
{
 if (useSTpackage  && !requireNamespace("SpatioTemporal", quietly = TRUE)) {
		useSTpackage <- FALSE
		warning("Argument 'useSTpackage' is TRUE, but SpatioTemporal package not available. Setting to FALSE.")
  }
  l.pars <- x
  m.mat <- R
  Sig <- block.Sig(l.pars,coords,reg.ind)
  beta <- inf.beta(l.pars, X= m.mat, coords, reg.ind, y, Sig=Sig, useSTpackage= useSTpackage)

  if(useSTpackage){
  	n.vec <- as.vector(by(reg.ind,reg.ind,length))
	Sig <- SpatioTemporal::makeCholBlock(Sig, n.blocks=length(n.vec), block.sizes=n.vec)
	Sig.inv <- SpatioTemporal:: invCholBlock(Sig, n.blocks=length(n.vec), block.sizes=n.vec)
  } else {
	Sig <- chol(Sig)
	Sig.inv <- chol2inv(Sig)
  }

 distval  <- stats::mahalanobis(x=y, center=m.mat%*%beta, cov=Sig.inv, inverted=T)
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


