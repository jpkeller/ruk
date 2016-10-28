###############################
##
## This file contains:
##
##	k.pred()
##	varcov.eff()
## 	block.Sig()
##	inf.beta()




## k.pred
##' @name k.pred
##' @title Kriging predictions
##
##' @description Internal function for computing kriging predictions (i.e. conditional expectations).
##
##' @details Not intended to be called directly. See \code{\link{predict.rlikfit}} for
##'		making kriging predictions from an \code{rlikfit} object.
## 
##' @param reg.dat matrix of regression data
##' @param dim.X TBA
##' @param beta TBA
##'	@param log.cov.pars TBA
##
##' @keywords internal
k.pred <- function(reg.dat,dim.X,beta,log.cov.pars)
{
  ind <- as.vector(by(reg.dat,reg.dat$p.ind,nrow))
  nm <- ind[1]
  nt <- sum(ind)
  y <- reg.dat$y.all[which(reg.dat$p.ind==1)]
  mu <- as.matrix(reg.dat[1:nm,3:(2+dim.X)])%*%beta
  vcov <- varcov.eff(reg.dat[,1:2],log.cov.pars)
  vcov12 <- vcov[1:nm,(nm+1):nt]
  vcov11 <- vcov[1:nm,1:nm]
  vcov22 <- vcov[(nm+1):nt,(nm+1):nt]
  iv11 <- chol(vcov11)
  iv11v12 <- backsolve(iv11, vcov12, transpose=TRUE)
  pred.var <- vcov22 - crossprod(iv11v12)
  iv11yM <- backsolve(iv11, y-mu, transpose=TRUE)
  reg.dat$y.all[which(reg.dat$p.ind==2)] <- as.matrix(reg.dat[(nm+1):nt,3:(2+dim.X)])%*%beta+crossprod(iv11v12, iv11yM)
#  pred.var <- vcov22-crossprod(vcov12,solve(vcov11,vcov12))
#  reg.dat$y.all[which(reg.dat$p.ind==2)] <- as.matrix(reg.dat[(nm+1):nt,3:(2+dim.X)])%*%beta+t(vcov12)%*%solve(vcov11)%*%(y-mu)
  reg.dat$v.all[which(reg.dat$p.ind==2)] <- diag(pred.var)
  reg.dat[,c('mp.ind','y.all','v.all')]
}



#================================
#CREATE EXPONENTIAL COVARIANCE MATRIX
#================================
##' @name varcov.eff
##' @title Creates exponential covariance matrix
##
##' @description Creates a covariance matrix according to an expontial model.
##
## details Details to be added.
##' @author Paul Sampson
##
##' @param coords Coordiantes
##' @param pars A vector of three log covariance parameters: nugget, sill, and range.
##
##' @return Returns an \eqn{n x n} matrix.
##' @importFrom stats dist
varcov.eff <- function(coords, pars)
{
  eye <- diag(1,nrow(coords))
  d <- as.matrix(stats::dist(coords))
  Sig <- exp(pars[1])*eye+exp(pars[2])*exp(-d/exp(pars[3]))
  Sig
}


##### TO do here:
# Remove the inversion of Sigma from this function,
# and streamline remaining code.

#=========================================
#Block diagonal spatial covariance matrix:
#=========================================
##' @name block.Sig
##' @title Creates block sigma
##
##' @description To do.
##
##' @details Details to be added.
##' @author Paul Sampson, Josh Keller
##
##' @param l.pars Log covariance parameters
##' @param coords Coordiantes
##' @param reg.ind Vector of region indicators
block.Sig <- function(l.pars, coords, reg.ind)
{
  coords <- coords[order(reg.ind),]
  reg.ind <- reg.ind[order(reg.ind)]
  n.vec <- as.vector(by(reg.ind,reg.ind,length))
  cum <- cumsum(n.vec)
  Sig <- matrix(0,nrow(coords),nrow(coords))
  Sig[1:n.vec[1],1:n.vec[1]]<- varcov.eff(coords[1:n.vec[1],],l.pars[1:3])
  if (length(n.vec)>1) {
    for (j in 2:length(cum)) {
        Sig[(cum[j-1]+1):(cum[j]),(cum[j-1]+1):(cum[j])] <- varcov.eff(coords[(cum[j-1]+1):(cum[j]),],l.pars[(3*j-2):(3*j)])
    }
  }
#  Sig.inv <- chol2inv(chol(Sig))
#  list('Sig.inv'=Sig.inv,'Sig'=Sig)
  Sig
}




#==============================
#FUNCTION TO INFER BETA FROM 
#MAXIMIZED COV PARS
#==============================
##' @name inf.beta
##' @title Infers Beta from maximized covariance pars.
##
##' @description Computes estimates of regression coefficients from
##'		covariance estimates.
##
##' @details This computes (X'SinvX)^(-1)(X'Sinv*y). Takes advantage of 
##'		functions in SpatioTemporal package if available.
##'
##' @author Paul Sampson, Casey Olives, Joshua Keller
##'
##' @param l.pars Log covariance parameters
##' @param X kriging covariates
##' @param coords Coordiantes
##' @param reg.ind Vector of region indicators
##' @param y outcome variable
##' @param Sig block-diagonal covariance matrix. Default value is \code{NULL}, and the
##'			matrix is computed by a call to \code{\link{block.Sig}}. 
##'			If it is already computed, providing it can speed up code.
##' @param useSTpackage See Details.
##
inf.beta<- function(l.pars, X, coords, reg.ind, y, Sig=NULL, useSTpackage=TRUE)
{
  if (useSTpackage  && !requireNamespace("SpatioTemporal", quietly = TRUE)) {
		useSTpackage <- FALSE
		warning("Argument 'useSTpackage' is TRUE, but SpatioTemporal package not available. Setting to FALSE.")
  }
  # to be sure, reorder by reg.ind as done elsewhere for block operations
  X <- X[order(reg.ind),, drop=FALSE]
  coords <- coords[order(reg.ind),]
  y <- y[order(reg.ind)]
  reg.ind <- reg.ind[order(reg.ind)]
  n.vec <- as.vector(by(reg.ind,reg.ind,length))
  if (is.null(Sig)){
  	Sig <- block.Sig(l.pars,coords,reg.ind)
  }
  if(useSTpackage){
	Sig <- SpatioTemporal::makeCholBlock(Sig, n.blocks=length(n.vec), block.sizes=n.vec)
	iSX <- SpatioTemporal::solveTriBlock(Sig, X, n.blocks=length(n.vec), block.sizes=n.vec, transpose=TRUE)
	iSy <- SpatioTemporal::solveTriBlock(Sig, y, n.blocks=length(n.vec), block.sizes=n.vec, transpose=TRUE)
	tXiSX <- crossprod(iSX)
	tXiSX <- chol(tXiSX)
	tXiSX.tXiS <- backsolve(tXiSX, t(iSX), transpose=TRUE)
	tXiSX.tXiS <- backsolve(tXiSX, tXiSX.tXiS, transpose=FALSE)
	beta <- tXiSX.tXiS %*% iSy
  } else {
	  Sig.inv <- chol2inv(chol(Sig))
	  beta <- solve(crossprod(X, Sig.inv %*% X), crossprod(X,Sig.inv %*% y))  	
  }
  beta
}

