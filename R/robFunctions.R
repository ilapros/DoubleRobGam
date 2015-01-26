

######## File written by Ilaria Prosdocimi (Ilaria.Prosdocimi@wis.kuleuven.be)



DevPois<- poisson()$dev.resid
DevBinom<-binomial()$dev.resid
DevGamma<-Gamma(link='log')$dev.resid  ### for now, only the log link
eexpp<-function(x) exp(x)/(1+exp(x))




###### everything that follows come straight from robustbase (also the comments - I left them, since they still hold)




### FIXME: should rather add these to R
pmin2 <- function(k,x) (x+k - abs(x-k))/2 
pmax2 <- function(k,x) (x+k + abs(x-k))/2



#### all the necessary functions for Epsi etc etc etc
## FIXME:  Do use a "robFamily", a  *list* of functions
## ------  which all have the same environment
##   ===> can get same efficiency as expressions, but better OO



### --- Normal -- family ---
EpsiNorm.init <- expression(
{
    E2f <- 0
})


EpsiNorm <- expression(
{
    0
})

Epsi2Norm <- expression(
{
    ## Calculation of E(psi^2) for the diagonal elements of A in matrix Q:
#     tcc^2*pnorm(mu-tcc*sqrt(phi),mean=mu,sd=sqrt(phi))+tcc^2*(1-pnorm(mu+tcc*sqrt(phi),mean=mu,sd=sqrt(phi)))+ (pnorm(mu+tcc*sqrt(phi),mean=mu,sd=sqrt(phi))-pnorm(mu-tcc*sqrt(phi),mean=mu,sd=sqrt(phi)))-tcc*exp(-0.5*tcc^2)*sqrt(2/pi)
tcc^2*pnorm(-tcc)+tcc^2*(1-pnorm(tcc))+(pnorm(tcc)-pnorm(-tcc))-tcc*exp(-0.5*tcc^2)*sqrt(2/pi)
})

EpsiSNorm <- expression(
{
    ## Calculation of E(psi*s) for the diagonal elements of B in the
    ## expression matrix M = 1/n t(X) %*% B %*% X:
   (pnorm(mu+tcc*sqrt(phi),mean=mu,sd=sqrt(phi))-pnorm(mu-tcc*sqrt(phi),mean=mu,sd=sqrt(phi)))/sqrt(phi)
})


Huber2.Norm <- function(phi,ns.resid,eta,Vmu,tcc,Epsi2val)
{
    sV <- sqrt(phi)
    snu <- 1/sqrt(phi)
#     eval(comp.Epsi.init)
#     compEpsi2 <- eval(Epsi2)
    sum(pmax2(-tcc,pmin(ns.resid*snu,tcc))^2)
}

DevNorm<-function (y, mu, wt)  wt * ((y - mu)^2)

phiNormEst.cl <- expression(
{
    ## Classical moment estimation of the dispersion parameter phi
    sum((y - mu)^2)/(nobs-ncoef)
})


### --- Poisson -- family ---

EpsiPois.init <- expression(
{
    dpH <- dpois(H, mu); dpH1 <- dpois(H-1, mu)
    dpK <- dpois(K, mu); dpK1 <- dpois(K-1, mu)
    pHm1 <- ppois(H-1, mu) ; pH <- pHm1 + dpH # = ppois(H,*)
    pKm1 <- ppois(K-1, mu) ; pK <- pKm1 + dpK # = ppois(K,*)
    E2f <- mu*(dpH1 - dpH - dpK1 + dpK) + pKm1 - pHm1
})

EpsiPois <- expression(
{
    tcc*(1 - pK - pH) + mu*(dpH - dpK)/sV
})

Epsi2Pois <- expression(
{
    ## Calculation of E(psi^2) for the diagonal elements of A in matrix Q:
    tcc^2 * (pH + 1 - pK) + E2f
})

EpsiSPois <- expression(
{
    ## Calculation of E(psi*s) for the diagonal elements of B in the
    ## expression matrix M = 1/n t(X) %*% B %*% X:
    tcc*(dpH + dpK) + E2f / sV
})

Huber2.Pois <- function(phi,ns.resid,eta,Vmu,tcc,Epsi2val)
{
    sV <- sqrt(Vmu*phi)
    snu <- 1/sqrt(phi)
#     eval(comp.Epsi.init)
#     compEpsi2 <- eval(Epsi2)
    sum(pmax2(-tcc,pmin(ns.resid*snu,tcc))^2) - sum(Epsi2val)
}



### --- Binomial -- family ---

EpsiBin.init <- expression({
    pK <- pbinom(K, ni, mu)
    pH <- pbinom(H, ni, mu)
    pKm1 <- pbinom(K-1, pmax2(0, ni-1), mu)
    pHm1 <- pbinom(H-1, pmax2(0, ni-1), mu)
    pKm2 <- pbinom(K-2, pmax2(0, ni-2), mu)
    pHm2 <- pbinom(H-2, pmax2(0, ni-2), mu)

    ## QlV = Q / V, where Q = Sum_j (j - mu_i)^2 * P[Y_i = j]
    ## i.e.  Q =	     Sum_j j(j-1)* P[.] +
    ##		 (1- 2*mu_i) Sum_j   j	 * P[.] +
    ##		     mu_i^2  Sum_j	   P[.]
    QlV <- mu/Vmu*(mu*ni*(pK-pH) +
		   (1 - 2*mu*ni) * ifelse(ni == 1, (H <= 0)*(K >= 1), pKm1 - pHm1) +
		   (ni - 1) * mu * ifelse(ni == 2, (H <= 1)*(K >= 2), pKm2 - pHm2))
})

EpsiBin <- expression(
{
    tcc*(1 - pK - pH) +
	ifelse(ni == 1, (- (H < 0) + (K >= 1) ) * sV,
	       (pKm1 - pHm1 - pK + pH) * mu * sni/sV)
})

Epsi2Bin <- expression(
{
    ## Calculation of E(psi^2) for the diagonal elements of A in matrix Q:
    tcc^2*(pH + 1 - pK) + QlV
})

EpsiSBin <- expression(
{
    ## Calculation of E(psi*s) for the diagonal elements of B in the
    ## expression matrix M = (X' B X)/n
    mu/Vmu*(tcc*(pH - ifelse(ni == 1, H >= 1, pHm1)) +
	    tcc*(pK - ifelse(ni == 1, K > 0,  pKm1))) + ifelse(ni == 0, 0, QlV / (sni*sV))
})

Huber2.Bin <- function(phi,ns.resid,eta,Vmu,tcc,Epsi2val)
{
    sV <- sqrt(Vmu*phi)
    snu <- 1/sqrt(phi)
#     eval(comp.Epsi.init)
#     compEpsi2 <- eval(Epsi2)
    sum(pmax2(-tcc,pmin(ns.resid*snu,tcc))^2) - sum(Epsi2val)
}


### --- Gamma -- family ---

Gmn <- function(t, nu) {
    ## Gm corrresponds to G * nu^((nu-1)/2) / Gamma(nu)
    snu <- sqrt(nu)
    snut <- snu+t
    r <- numeric(length(snut))
    ok <- snut > 0
    r[ok] <- {
	nu <- nu[ok]; snu <- snu[ok]; snut <- snut[ok]
	exp((nu-1)/2*log(nu) - lgamma(nu) - snu*snut + nu*log(snut))
    }
    r
}

EpsiGamma.init <- expression({
    nu <- 1/phi      ## form parameter nu
    snu <- 1/sqrt(phi) ## sqrt (nu)
    pPtc <- pgamma(snu + c(-tcc,tcc), shape=nu, rate=snu)
    pMtc <- pPtc[1]
    pPtc <- pPtc[2]
    aux2 <- tcc*snu
    GLtcc <- Gmn(-tcc,nu)
    GUtcc <- Gmn( tcc,nu)
})

EpsiGamma <- expression( tcc*(1-pPtc-pMtc) + GLtcc - GUtcc )

EpsiSGamma <- expression( ((GLtcc - GUtcc) + snu*(pPtc-pMtc))/mu )

Epsi2Gamma <- expression({
    (tcc^2*(pMtc+1-pPtc)+ (pPtc-pMtc) +
     (GLtcc*(1-aux2) - GUtcc*(1+aux2))/snu )
})


phiGammaEst.cl <- expression(
{
    ## Classical moment estimation of the dispersion parameter phi
    sum(((y - mu)/mu)^2)/(nobs-ncoef)
})

phiGammaEst <- expression(
{
    ## robust estimation of the dispersion parameter by
    ## Huber's porposal 2
    sphi <- uniroot(Huberprop2, interval=range(residP^2),
                    ns.resid=residP, mu=mu, Vmu=Vmu, tcc=tcc)$root
})

if(FALSE) ## ...  (MM ?? FIXME : use  EpsiGamma.init(), Epsi2Gamma() !! )
Huberprop2.gehtnicht <- function(phi, ns.resid, mu, Vmu, tcc)
{
    sV <- sqrt(Vmu*phi)
    H <- floor(mu - tcc* sV)
    K <- floor(mu + tcc* sV)
    nobs <- length(mu)
    ##
    eval(Epsi.init)
    compEpsi2 <- eval(Epsi2)
    ##
    ## return h :=
    sum(pmax2(-tcc,pmin(ns.resid*snu,tcc))^2) -  nobs*compEpsi2
}

Huberprop2 <- function(phi,ns.resid,mu,Vmu,tcc)
{
    sV <- sqrt(Vmu*phi)
    H <- floor(mu - tcc* sV)
    K <- floor(mu + tcc* sV)
    nobs <- length(mu)

    nu <- 1/phi         ## form parameter  nu
    snu <- 1/sqrt(phi)  ## sqrt (nu)
    pPtc <- pgamma(snu + c(-tcc,tcc), shape=nu, rate=snu)
    pMtc <- pPtc[1]
    pPtc <- pPtc[2]

    ts <- tcc*snu
    GLtcc <- Gmn(-tcc,nu) *(1-ts)/snu
    GUtcc <- Gmn( tcc,nu) *(1+ts)/snu
    ##
    compEpsi2 <- tcc^2 + (pPtc - pMtc)*(1-tcc^2) + GLtcc - GUtcc
    ## return h :=
    sum(pmax2(-tcc,pmin(ns.resid*snu,tcc))^2) -  nobs*compEpsi2
}






######## weights for leverage points

wts_HiiDist <- function(X) {
    x <- qr(X)
    Hii <- rowSums(qr.qy(x, diag(1, nrow = NROW(X), ncol = x$rank))^2)
    sqrt(1-Hii)
}

wts_HiiDist <- function(X) {
    x <- qr(X)
    Hii <- rowSums(qr.qy(x, diag(1, nrow = NROW(X), ncol = x$rank))^2)
    s<-sqrt(1-Hii);s[is.na(s)]<-0
    s
}



wts_RobDist <- function(X, intercept, covFun)
{
    if(intercept) {
	X <- as.matrix(X[, -1])
	Xrc <- covFun(X)
	dist2 <- mahalanobis(X, center = Xrc$center, cov = Xrc$cov)
    }
    else {
	if(!is.matrix(X)) X <- as.matrix(X)
	Xrc <- covFun(X)
	mu <- as.matrix(Xrc$center)
	Mu <- Xrc$cov + tcrossprod(mu)
	dist2 <- mahalanobis(X, center = rep(0,ncol(X)), cov = Mu)
    }
    ncoef <- ncol(X) ## E[chi^2_p] = p
    1/sqrt(1+ pmax2(0, 8*(dist2 - ncoef)/sqrt(2*ncoef)))
}


