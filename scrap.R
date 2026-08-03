TODO

1. x to y
2. move supp.x into params in ks.study
3. rename GH$x to GH$theta.grid

## DEBUGGING

## check variance of mle's to see if crazy k.vals could be right. no, variances are about .2, something is wrong with k formulas


## set.seed(1)
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
## n <- 10
## ns <- 5:15
## by.n <- sapply(ns, function(n) {
## by.n <- parallel::mclapply(ns,mc.cores=4,FUN= function(n) {
## cat('.')
m.T <- rep(20,n)
m.C <- rep(15,n)
s <- sapply(1:n,function(i)sample(max(m.T[i],m.C[i]),1))
## s <- rep(10,n)
n <- 5
stats <- replicate(round(n^2.3), {
    theta <- rnorm(n)
    odds <- exp(mu.true+sqrt(tau2.true)*theta)
    x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
    params <- cbind(m.T=m.T,m.C=m.C,s=s)
    LR <- LR.gaussian(x,params,mu.null,loglik.unit.hyperg)$mle.full
})


source('utils.R')
mma.str <- inf.observed.str <- '(E^(x*(theta*Sqrt[beta[1]] + beta[2]))*Binomial[mC, s - x]*Binomial[mT, x]*
  (x^2*Inactive[Sum][E^(j*(theta*Sqrt[beta[1]] + beta[2]))*Binomial[mC, -j + s]*
       Binomial[mT, j], {j, 0, mT}]^2 + 
   2*Inactive[Sum][E^(j*(theta*Sqrt[beta[1]] + beta[2]))*j*Binomial[mC, -j + s]*
       Binomial[mT, j], {j, 0, mT}]^2 - 
   Inactive[Sum][E^(j*(theta*Sqrt[beta[1]] + beta[2]))*Binomial[mC, -j + s]*Binomial[mT, j], 
     {j, 0, mT}]*(2*x*Inactive[Sum][E^(j*(theta*Sqrt[beta[1]] + beta[2]))*j*
        Binomial[mC, -j + s]*Binomial[mT, j], {j, 0, mT}] + 
     Inactive[Sum][E^(j*(theta*Sqrt[beta[1]] + beta[2]))*j^2*Binomial[mC, -j + s]*
       Binomial[mT, j], {j, 0, mT}])))/
 Inactive[Sum][E^(j*(theta*Sqrt[beta[1]] + beta[2]))*Binomial[mC, -j + s]*Binomial[mT, j], 
   {j, 0, mT}]^3'
mma.str <- gsub('Inactive[Sum]','Sum',mma.str,fixed=TRUE)
mma.str <- gsub('mC','m.C',mma.str,fixed=TRUE)
mma.str <- gsub('mT','m.T',mma.str,fixed=TRUE)
mma.str <-  mma.to.R(mma.str,
                     funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","Hypergeometric2F1"),
                     rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",Hypergeometric2F1="gsl::hyperg_2F1"))
finf <- 
    str.to.fun(mma.str,
               args = c('beta','theta','x','m.T','m.C','s'),
               binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
x <- 5;m <- m.T <- 20; m.C <- 15; s <- 10;
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
beta <- c(mu.null,tau2.null)
GH <- fastGHQuad::gaussHermiteData(30)
GH$x <- sqrt(2)*GH$x
GH$w <- GH$w / sqrt(pi)
names(GH) <- c('theta.grid','w')
finf2 <- function(beta,x,m.T,m.C,s) sum(GH$w*sapply(GH$theta.grid,function(theta)finf(beta,theta,x,m.T,m.C,s)))
## finf2(beta,x,m.T,m.C,s)
pmf <- function(theta,x) { # numerically unsafe
    numer <- choose(m.T,x)*choose(m.C,s-x)*exp((mu.null+sqrt(tau2.null)*theta)*x)
    supp <- max(0,s-m.C):min(s,m.T)
    denom <- rowSums(  outer(theta,supp,function(thetaj,xj)choose(m.T,xj)*choose(m.C,s-xj)*exp((mu.null+sqrt(tau2.null)*thetaj)*xj) ) )
    ## denom <- rowSums(  outer(theta,0:m.Ti,function(thetaj,xj)choose(m.Ti,xj)*choose(m.Ci,si-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
    numer / denom
}
pmf2 <- function(x)sum(GH$w*pmf(GH$theta.grid,x))
supp <- max(0,s-m.C):min(s,m.T)
sum(sapply(supp,function(x)finf2(beta,x,m.T,m.C,s)*pmf2(x)))






## poisson-normal model



## debugging--variances implied by k values seem off


## expected score is 0 -- fine
source('utils.R')
mma.str <- '-(((E^(T*(theta*Sqrt[beta[1]] + beta[2])))^x*T*(E^(T*(theta*Sqrt[beta[1]] + beta[2])) - x))/
  (E^E^(T*(theta*Sqrt[beta[1]] + beta[2]))*factorial[x]))'
mma.str <-  mma.to.R(mma.str,
                     funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial"),
                     rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
fun <- 
    str.to.fun(mma.str,
               args = c('beta','theta','x','T'),
               binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
GH <- fastGHQuad::gaussHermiteData(100)
GH$x <- sqrt(2)*GH$x
GH$w <- GH$w / sqrt(pi)
names(GH) <- c('theta.grid','w')
fun2 <- function(theta,x)fun(beta=c(1/2,0),theta=theta,x=x,T=1)
mat <- outer(GH$theta.grid,0:50,fun2)
sum(GH$w*mat)

## information for mu and tau2
source('utils.R')
mma.str.0 <- '(E^(theta*Sqrt[beta[1]] + beta[2])*T)^x/(E^(E^(theta*Sqrt[beta[1]] + beta[2])*T)*factorial[x])'
## mma.str.1 <- '((E^(theta*Sqrt[beta[1]] + beta[2])*T)^x*(-(E^(theta*Sqrt[beta[1]] + beta[2])*T) + x))/
##  (E^(E^(theta*Sqrt[beta[1]] + beta[2])*T)*factorial[x])'
## mma.str.2 <- '((E^(theta*Sqrt[beta[1]] + beta[2])*T)^x*(E^(2*(theta*Sqrt[beta[1]] + beta[2]))*T^2 + x^2 - 
##    E^(theta*Sqrt[beta[1]] + beta[2])*T*(1 + 2*x)))/(E^(E^(theta*Sqrt[beta[1]] + beta[2])*T)*
##   factorial[x])'
mma.str.1 <- '-1/2*((E^(theta*Sqrt[beta[1]] + beta[2])*T)^x*theta*(E^(theta*Sqrt[beta[1]] + beta[2])*T - 
    x))/(E^(E^(theta*Sqrt[beta[1]] + beta[2])*T)*Sqrt[beta[1]]*factorial[x])'
mma.str.2 <- '((E^(theta*Sqrt[beta[1]] + beta[2])*T)^x*theta*(x*(-1 + theta*x*Sqrt[beta[1]]) - 
   E^(theta*Sqrt[beta[1]] + beta[2])*T*(-1 + theta*(1 + 2*x)*Sqrt[beta[1]]) + 
   E^(2*(theta*Sqrt[beta[1]] + beta[2]))*T^2*theta*Sqrt[beta[1]]))/
 (4*E^(E^(theta*Sqrt[beta[1]] + beta[2])*T)*beta[1]^(3/2)*factorial[x])'
fsums <- sapply(c('0'=mma.str.0,'1'=mma.str.1,'2'=mma.str.2), function(mma.str) {
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
    fun <- 
        str.to.fun(mma.str,
                   args = c('beta','theta','x','T'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
    GH <- fastGHQuad::gaussHermiteData(100)
    GH$x <- sqrt(2)*GH$x
    GH$w <- GH$w / sqrt(pi)
    names(GH) <- c('theta.grid','w')
    fun2 <- function(theta,x)fun(beta=c(1/2,0),theta=theta,x=x,T=1)
    mat <- outer(GH$theta.grid,0:50,fun2)
    ## sum(GH$w*mat)
    ## GH$w%*%mat
    colSums(GH$w*mat)
})
sum((fsums[,'2'] * fsums[,'0'] - fsums[,'1']^2) / fsums[,'0'])
## > [1] -0.7201051
##  [1] -0.3897629

## off diagonal
source('utils.R')
mma.str.0 <- '(E^(theta*Sqrt[beta[1]] + beta[2])*T)^x/(E^(E^(theta*Sqrt[beta[1]] + beta[2])*T)*factorial[x])'
mma.str.1 <- '-1/2*((E^(theta*Sqrt[beta[1]] + beta[2])*T)^x*theta*(E^(theta*Sqrt[beta[1]] + beta[2])*T - 
    x))/(E^(E^(theta*Sqrt[beta[1]] + beta[2])*T)*Sqrt[beta[1]]*factorial[x])'
mma.str.2 <- '((E^(theta*Sqrt[beta[1]] + beta[2])*T)^x*(-(E^(theta*Sqrt[beta[1]] + beta[2])*T) + x))/
 (E^(E^(theta*Sqrt[beta[1]] + beta[2])*T)*factorial[x])'
mma.str.12 <- '((E^(theta*Sqrt[beta[1]] + beta[2])*T)^x*theta*(E^(2*(theta*Sqrt[beta[1]] + beta[2]))*T^2 + 
   x^2 - E^(theta*Sqrt[beta[1]] + beta[2])*T*(1 + 2*x)))/
 (2*E^(E^(theta*Sqrt[beta[1]] + beta[2])*T)*Sqrt[beta[1]]*factorial[x])'
    fsums <- sapply(c('0'=mma.str.0,'1'=mma.str.1,'2'=mma.str.2,'12'=mma.str.12), function(mma.str) {
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial"))
    fun <- 
        str.to.fun(mma.str,
                   args = c('beta','theta','x','T'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
    GH <- fastGHQuad::gaussHermiteData(100)
    GH$x <- sqrt(2)*GH$x
    GH$w <- GH$w / sqrt(pi)
    names(GH) <- c('theta.grid','w')
    fun2 <- function(theta,x)fun(beta=c(1/2,0),theta=theta,x=x,T=1)
    mat <- outer(GH$theta.grid,0:50,fun2)
    ## sum(GH$w*mat)
    ## GH$w%*%mat
    colSums(GH$w*mat)
})
sum((fsums[,'12'] * fsums[,'0'] - fsums[,'1']*fsums[,'2']) / fsums[,'0'])
## > [1] -0.7201051
##  [1] -0.3897629





## REFACTORED CODE


loglik.unit.hyperg <- function(mu,tau2,xi,paramsi,GH=NULL) {
    m.Ti <- paramsi[,'m.T']
    m.Ci <- paramsi[,'m.C']
    si <- paramsi[,'s']
    ## pmf <- Vectorize(function(theta)BiasedUrn::dFNCHypergeo(xi,m.Ti,m.Ci,si,exp(mu+sqrt(tau2)*theta))) # slow
    pmf <- function(theta) { # numerically unsafe
        numer <- choose(m.Ti,xi)*choose(m.Ci,si-xi)*exp((mu+sqrt(tau2)*theta)*xi)
        supp <- max(0,si-m.Ci):min(si,m.Ti)
        denom <- rowSums(  outer(theta,supp,function(thetaj,xj)choose(m.Ti,xj)*choose(m.Ci,si-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
        ## denom <- rowSums(  outer(theta,0:m.Ti,function(thetaj,xj)choose(m.Ti,xj)*choose(m.Ci,si-xj)*exp((mu+sqrt(tau2)*thetaj)*xj) ) )
        numer / denom
    }
    ## pmf <- function(theta) {
    ##     ## logpmf.raw <- function(xi,theta)lgamma(m.Ti+1) - lgamma(xi+1) -lgamma(m.Ti-xi+1) + lgamma(m.Ci+1) - lgamma(si-xi+1) - lgamma(m.Ci-si+xi+1) + (mu+sqrt(tau2)*theta)*xi
    ##     ## lognumer <- logpmf.raw(xi,theta)
    ##     ## logpmf.mat <- outer(0:m.Ti,theta,logpmf.raw) # speed this up
    ##     logpmf.part <- function(xi)lgamma(m.Ti+1) - lgamma(xi+1) -lgamma(m.Ti-xi+1) + lgamma(m.Ci+1) - lgamma(si-xi+1) - lgamma(m.Ci-si+xi+1)
    ##     lognumer <- logpmf.part(xi)+ (mu+sqrt(tau2)*theta)*xi
    ##     ## logpmf.mat <- logpmf.part(0:m.Ti) + outer(0:m.Ti,theta)*sqrt(tau2) + (0:m.Ti)*mu
    ##     supp <- max(0,si-m.Ci):min(si,m.Ti)
    ##     logpmf.mat <- logpmf.part(supp) + outer(supp,theta)*sqrt(tau2) + supp*mu
    ##     colmaxes <- apply(logpmf.mat,2,max)        
    ##     logdenom <- colmaxes + log(rowSums(exp(t(logpmf.mat) - colmaxes)))
    ##     exp(lognumer - logdenom)
    ## }
    integral.approx <- if(is.null(GH)) integrate(function(theta)pmf(theta)*dnorm(theta),-Inf,Inf)$val
    else
        sum(GH$w*pmf(GH$x))
    return(log(integral.approx))
}



loglik.unit <- function(mu,tau2,xi,paramsi,cond.pmf,GH=NULL) {
    ## browser()
    if(is.null(GH)) {
        likelihood <- integrate(function(theta)cond.pmf(xi,theta,mu,tau2,paramsi)*dnorm(theta),-Inf,Inf)$val 
    } else {
        likelihood <- sum(GH$w*cond.pmf(xi,GH$x,mu,tau2,paramsi))
    }
    return(log(likelihood))
}
## set.seed(1)
## mu <- runif(1); tau2 <- rexp(1); paramsi <- c(T=10); xi <- sample(paramsi['T'],1)
## GH <- fastGHQuad::gaussHermiteData(30)
## GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
## loglik.unit(mu,tau2,xi,paramsi,cond.pmf.poisson,GH=GH)




## ks.study <- function(m,pmf.derivs.fun.ma,ks.fun,GH=NULL) {
##     if(is.null(GH)) {
##         GH <- fastGHQuad::gaussHermiteData(40)
##         GH$x <- sqrt(2)*GH$x
##         GH$w <- GH$w / sqrt(pi)
##         names(GH) <- c('theta.grid','w')
##     }
##     ## sum ma-level g's over theta, getting a vector g(x,m) with
##     ## unit-level components
##     pmf.derivs.unit <- lapply(pmf.derivs.fun.ma, function(pmf.deriv.unit) {
##         ## pmf.deriv.unit <- Vectorize(pmf.deriv.unit,c('theta','x'))
##         supp <- 0:50 #!!!!!!!!
##         g.mat <- outer(X=GH$theta.grid,Y=supp,function(theta,x)pmf.deriv.unit(theta,x,T)) # these are globals!!!
##         g.unit <- colSums(GH$w * g.mat)
##     })
##     ## apply ks.fun to unit-level g's to get unit-level k's, take sum
##     ## weighted by binomial coefficient to get study-level k
##     k.study <- sapply(seq_along(ks.fun), function(i) {
##         k.fun <- ks.fun[[i]]$fun; args <- ks.fun[[i]]$args
##         k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
##         k.val <- sum( k.fun.unit, na.rm=TRUE) #!!!!
##     })
##     return(k.study)
## }



## ## read in pmf derivatives
## source('utils.R')
## mma.str.path <- './mma-formulas/pmf-hyperg/'
## filelist <- dir(mma.str.path)
## filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## ## file <- './mma-formulas/quad/g/g-10.txt'
## pmf.derivs.fun <- lapply(filelist, function(file) {
##     mma.str <-paste(readLines(file,warn=FALSE),collapse='')
##     ## mma.str <- gsub('Pi','pi',mma.str)
##     mma.str <- gsub('Inactive[Sum]','Sum',mma.str,fixed=TRUE)
##     mma.str <- gsub('mC','m.C',mma.str,fixed=TRUE)
##     mma.str <- gsub('mT','m.T',mma.str,fixed=TRUE)
##     mma.str <-  mma.to.R(mma.str,
##                          funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","Hypergeometric2F1"),
##                          rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",Hypergeometric2F1="gsl::hyperg_2F1"))
##     g.deriv <- 
##         str.to.fun(mma.str,
##                    args = c('beta','theta','x','m.T','m.C','s'),
##                    binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt))
## })
## pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
## names(pmf.derivs.fun) <- pmf.deriv.orders
## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,m.T,m.C,s)pmf.deriv.fun(beta,theta,x,m.T,m.C,s))
## ## need to vectorize these generically
## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun.ma,function(pmf.deriv.fun.ma) Vectorize(pmf.deriv.fun.ma,c('theta','x')))





DEBUG



## checking poisson model numerical issues
## 8c-4 poisson
source('utils.R')
mu.null <- 0
tau2.null <- 1/2
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
T <- 50
supp.x <- 0:1e2
params <- c(T=T)
## ks.df <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params,supp.x=supp.x,mu.null=mu.null,tau2.null=tau2.null)
## bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
## dput(bartlett.C)
beta <- c(tau2.null,mu.null)
## if(is.null(GH)) {
    GH <- fastGHQuad::gaussHermiteData(50)
    GH$x <- sqrt(2)*GH$x
    GH$w <- GH$w / sqrt(pi)
    names(GH) <- c('theta.grid','w')
## }

pmf.deriv.fun <- pmf.derivs.fun[[1]]
pmf.derivs.unit <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
    ## supp <- max(0,s-m.C):min(s,m.T)
    ## g.mat <- outer(X=GH$theta.grid,Y=supp,function(theta,x)pmf.deriv.fun.ma(theta,x,m,m.C,s)) # these are globals!!!
    g.mat <- outer(X=GH$theta.grid,Y=supp.x,function(theta,x)pmf.deriv.fun(beta=beta,theta=theta,x=x,params=params)) 
    g.unit <- colSums(GH$w * g.mat)
})
max(sapply(pmf.derivs.unit,function(v)max(v)))
min(sapply(pmf.derivs.unit,function(v)min(v)))


(E^(theta * sqrt(beta[1]) + beta[2]) * params["T"])^x/(E^(E^(theta * 
    sqrt(beta[1]) + beta[2]) * params["T"]) * factorial(x))




## read in pmf derivatives
source('utils.R')
mu.null <- 0
tau2.null <- 1/2
ks.fun <- readRDS('ks-fun.rds')

mma.str.path <- './mma-formulas/pmf-poisson/'
filelist <- dir(mma.str.path)
filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
## file <- './mma-formulas/quad/g/g-10.txt'
pmf.derivs.fun <- lapply(filelist, function(file) {
    mma.str <-paste(readLines(file,warn=FALSE),collapse='')
    ## mma.str <- gsub('Pi','pi',mma.str)
    ## mma.str <- gsub('Inactive[Sum]','Sum',mma.str,fixed=TRUE)
    ## mma.str <- gsub('mC','m.C',mma.str,fixed=TRUE)
    ## mma.str <- gsub('mT','m.T',mma.str,fixed=TRUE)
    mma.str <-  mma.to.R(mma.str,
                         funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial","factorial","dpois"),
                         rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose",factorial="factorial",dpois="dpois"))
    g.deriv <- 
        str.to.fun(mma.str,
                   args = c('beta','theta','x','params'),
                   binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt,dpois=dpois))
})
pmf.deriv.orders <- paste0('pmfSum',substr(filelist,nchar(filelist)-5,nchar(filelist)-4))
names(pmf.derivs.fun) <- pmf.deriv.orders
## saveRDS(pmf.derivs.fun,file='pmf-derivs-fun-poisson.rds',compress=FALSE)
## pmf.derivs.fun.ma <- lapply(pmf.derivs.fun,function(pmf.deriv.fun)function(theta,x,params)pmf.deriv.fun(beta,theta,x,params))



## pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
T <- 1
supp.x <- 0:2e2
params <- c(T=T)

ks.df <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params,supp.x=supp.x,mu.null=mu.null,tau2.null=tau2.null)
bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
dput(bartlett.C)




## check z-stats with new poisson using recurrence 
source('utils.R')
mu.null <- -4
tau2.null <- .2
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
T <- 1e3
supp.x <- 0:1e3
params <- c(T=T)
ks.df <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params,supp.x=supp.x,mu.null=mu.null,tau2.null=tau2.null,n.GH=1e2)
k.22 <- ks.df[ks.df[,'inner']=='22'  & ks.df[,'outer']=='0','k.val']
k.11 <- ks.df[ks.df[,'inner']=='11'  & ks.df[,'outer']=='0','k.val']
k.12 <- ks.df[ks.df[,'inner']=='12'  & ks.df[,'outer']=='0','k.val']
print(k.22); print(k.11)
bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
dput(bartlett.C)


n <- 5e1
inf <- -matrix(c(k.11,k.12,k.12,k.22),nrow=2)
vcov.hat <- solve(inf)
## T <- rep(1,n)
z.stats <- replicate(1e2, {
    theta <- rnorm(n)
    ## T <- 4 + rexp(n)
    lambda <- T*exp(mu.null+sqrt(tau2.null)*theta)
    x <- rpois(n,lambda)
    params <- cbind(T=rep(T,n))
    mles <- LR.gaussian(x,params,mu.null,cond.pmf.poisson,n.GH=100)$mle.full
    (rev(mles)-c(tau2.null,mu.null)) / sqrt(diag(vcov.hat)/n)
})
op <- par(mfrow=c(1,2))
qqnorm(z.stats['mu',]); abline(0,1)
qqnorm(z.stats['tau2',]); abline(0,1)
par(op)



n <- 5
stats <- replicate(5e2, {
    theta <- rnorm(n)
    ## T <- 4 + rexp(n)
    lambda <- T*exp(mu.null+sqrt(tau2.null)*theta)
    x <- rpois(n,lambda)
    params <- cbind(T=rep(T,n))
    LR <- LR.gaussian(x,params,mu.null,cond.pmf.poisson,n.GH=1e2)$LR
    ## LR <- LR.gaussian(x,params,mu.null,loglik.unit.hyperg)$LR
    LR.bartlett <- LR / (1+bartlett.C/n)
    c(LR=LR,LR.bartlett=LR.bartlett)
})
rowMeans(stats-1)
## hist(stats['LR',],xlim=c(0,2))
## abline(v=1,col=2)
## abline(v=mean(stats['LR',]),col=3)
## abline(v=mean(stats['LR.bartlett',]),col=3,lty=2)



DATA ENTRY

## niel.weise <- data.frame(
##   study = 1:18,
##   Y0    = c(3, 3, 9, 7, 6, 4, 3, 2, 19, 2, 7, 1, 5, 11, 0, 1, 3, 1),
##   n0    = c(117, 35, 195, 136, 157, 139, 177, 39, 103, 122, 64, 58, 175, 180, 105, 262, 362, 69),
##   Y1    = c(0, 1, 2, 0, 5, 1, 1, 1, 1, 1, 0, 0, 3, 6, 0, 0, 1, 4),
##   n1    = c(116, 44, 208, 130, 151, 98, 174, 74, 97, 113, 66, 70, 188, 187, 118, 252, 345, 64)
## )
## write.table(niel.weise,'niel-weise.txt')

## mu.null <- with(niel.weise,qlogis(mean(Y0/n0)))
## source('utils.R')
## LR.gaussian(x=niel.weise$Y0,params=cbind(m=niel.weise$n0),mu.null=0,cond.pmf=cond.pmf.binomial.prop,n.GH=50)














OLD VERSIONS OF UTILS


## ## k estimates for a single study
## ## pmf.deriv.fun.ma should contain functions vectorized in theta and
## ## x.
## ks.study <- function(pmf.derivs.fun.ma,ks.fun,params,support,GH=NULL) {
##     ## browser()
##     if(is.null(GH)) {
##         GH <- fastGHQuad::gaussHermiteData(50)
##         GH$x <- sqrt(2)*GH$x
##         GH$w <- GH$w / sqrt(pi)
##         names(GH) <- c('theta.grid','w')
##     }
##     ## sum ma-level g's over theta, getting a vector g(xi) with
##     ## unit-level components. 
##     pmf.derivs.unit <- lapply(pmf.derivs.fun.ma, function(pmf.deriv.fun.ma) {
##         ## supp <- max(0,s-m.C):min(s,m.T)
##         ## g.mat <- outer(X=GH$theta.grid,Y=supp,function(theta,x)pmf.deriv.fun.ma(theta,x,m,m.C,s)) # these are globals!!!
##         g.mat <- outer(X=GH$theta.grid,Y=support,function(theta,x)pmf.deriv.fun.ma(theta,x,params)) # these are globals!!!
##         g.unit <- colSums(GH$w * g.mat)
##     })
##     ## apply ks.fun to unit-level g's to get unit-level k's, take sum
##     ## weighted by binomial coefficient to get study-level k
##     k.study <- sapply(seq_along(ks.fun), function(i) {
##         k.fun <- ks.fun[[i]]$fun; args <- ks.fun[[i]]$args
##         k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
##         k.val <- sum( k.fun.unit)
##     })
##     return(k.study)
## }


## ## allowing for theta-only parts of the pmf derivative to be passed
## ## in. rework ks.fun to store order and args as attributes not list
## ## elements.
## ## update: this approach was abandoned in favor of better formulas (using recurrences)
## ks.study <- function(pmf.derivs.fun,ks.fun,theta.parts.fun=list(),params,supp.x,mu.null,tau2.null,GH=NULL) {
##     ## browser()
##     beta <- c(tau2.null,mu.null)
##     if(is.null(GH)) {
##         GH <- fastGHQuad::gaussHermiteData(50)
##         GH$x <- sqrt(2)*GH$x
##         GH$w <- GH$w / sqrt(pi)
##         names(GH) <- c('theta.grid','w')
##     }
##     theta.parts <- lapply(theta.parts.fun, function(f) {
##         theta.part.mat <- outer(GH$theta.grid,supp.x,function(theta,x)f(beta,theta,x,params))
##         rowSums(theta.part.mat)
##     })
##     pmf.derivs.unit <- lapply(pmf.derivs.fun, function(f) {
##         arg.list <- c(list(
##             theta = rep(GH$theta.grid,length(supp.x)),
##             x = rep(supp.x,each=length(GH$theta.grid))),
##             lapply(theta.parts[attr(f,'args')], function(v)rep(v,length(supp.x))),
##             beta=list(beta),
##             params=list(params)
##             )
##         pmf.deriv.unit <- do.call(f,arg.list)
##         dim(pmf.deriv.unit) <- c(length(GH$theta.grid),length(supp.x))
##         GH$w %*% pmf.deriv.unit
##     })
##     k.vals <- sapply(seq_along(ks.fun), function(i) {
##         k.fun <- ks.fun[[i]]
##         args <- attr(k.fun,'args')
##         k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
##         k.val <- sum( k.fun.unit )
##     })
##     k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
##     ## k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
##     ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
##     return(ks.df)
## }




## LR.gaussian <- function(x,params,mu.null,loglik.unit,n.GH=30) {
##         GH <- fastGHQuad::gaussHermiteData(30)
##         GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
##         loglik <- function(mu,tau2,x,params,GH=NULL) sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,,drop=FALSE],GH=GH)))
##         logtau2.mle.redu <- optimize(function(logtau2)loglik(mu.null,exp(logtau2),x,params,GH),interval=c(-5,5),maximum=TRUE)$max
##         tau2.mle.redu <- exp(logtau2.mle.redu)
##         mle.full <- optim(par=c(mu.null,log(tau2.mle.redu)),fn=function(par)-loglik(par[1],exp(par[2]),x,params,GH),method='BFGS')$par
##         mu.mle.full <- mle.full[1]; tau2.mle.full <- exp(mle.full[2])
##         LR <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.null,tau2.mle.redu,x,params,GH))
##         return(list(mle.full=c(mu=mu.mle.full,tau2=tau2.mle.full),mle.redu=c(mu=mu.null,tau2=tau2.mle.redu), LR=LR))
## }



read.mma <- function(mma.str.path) {
    filelist <- dir(mma.str.path)
    filelist <- paste0(mma.str.path,filelist[grep('\\.txt',filelist)])
    mma.strings <- sapply(filelist,function(file)paste(readLines(file,warn=FALSE),collapse=''))
    mma.strings <- sapply(mma.strings,
                          function(mma.str) mma.to.R(mma.str,
                                                     funs = c("Sqrt", "Log", "LogisticSigmoid", "Binomial"),
                                                     rename = c(Sqrt = "sqrt",Log = "log",LogisticSigmoid = "plogis",Binomial = "choose")
                                                     )
                          )
    l.derivs <- sapply(mma.strings,function(mma.str) 
        str.to.fun(
            gsub(pattern=c('m\\[j\\]'),replacement=c('m'),mma.str),
            args = c('theta','m','x','beta'),
            binding=list(E=exp(1),plogis=plogis,log=log,sqrt=sqrt)
        ))
}


## bartlett.ma <- function(params,ks.fun,pmf.derivs.fun) {   
##     ks <- lapply(1:nrow(params), function(i) {
##         supp.x <- params[i,'x.min']:params[i,'x.max']
##         ks.df.study <- ks.study0(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params[i,],mu.null=mu.null,tau2.null=tau2.null) 
##     })
##     k.vals.ma <- sapply(ks,function(df)df$k.val) |> rowSums()
##     ks.df.ma <- data.frame(inner=ks[[1]]$inner, outer=ks[[1]]$outer,k.val=k.vals.ma)
##     bartlett.C <- lawley(ks.df.ma,2) - lawley(ks.df.ma,1)
##     return(bartlett.C)
## }
## params <- cbind(m=ms,x.min=0,x.max=ms)
## bartlett.ma(params,ks.fun,pmf.derivs.fun)


DATA INPUT

stijnen.vi <- data.frame(
  Study = 1:9,
  Y1    = c(7, 8, 3, 6, 1, 1, 17, 3, 2),
  T1    = c(1491, 1638, 11484, 1446, 370, 786, 6840, 1156, 400),
  Y0    = c(11, 8, 14, 10, 1, 8, 15, 7, 3),
  T0    = c(1988, 1460, 10962, 1502, 482, 913, 6840, 1015, 440)
)
write.table(stijnen.vi,file='stijnen-table-VI.txt')


stijnen.v <- data.frame(
  Study = 1:9,
  TP    = c(33, 18, 20, 19, 44, 18, 27, 21, 18), 
  FN    = c(6, 6, 13, 0, 6, 3, 1, 0, 1),        
  n1    = c(39, 24, 33, 19, 50, 21, 28, 21, 19),
  FP    = c(5, 5, 0, 0, 10, 1, 4, 1, 1),        
  TN    = c(35, 10, 41, 19, 19, 9, 21, 9, 20),  
  n0    = c(40, 15, 41, 19, 29, 10, 25, 10, 21) 
)
write.table(stijnen.v,'stijnen-table-V.txt')

# Set common aliases matching metafor / Stijnen notation:
patwardhan$Y1 <- patwardhan$TP
patwardhan$Y0 <- patwardhan$FP


DEBUGGING POISSON PROFILE CURVE. TURNED OUT TO BE SUPPORT WAS TOO SMALL.



## DEBUG poisson profile curve

mu.profile.grid <- c(-0.00568896507345706, 0.142381275751646, 0.29045151657675, 
0.438521757401854, 0.586591998226957, 0.734662239052061, 0.882732479877164, 
1.03080272070227, 1.17887296152737, 1.32694320235247, 1.47501344317758, 
1.62308368400268, 1.77115392482779, 1.91922416565289, 2.06729440647799, 
2.2153646473031, 2.3634348881282, 2.5115051289533, 2.65957536977841, 
2.80764561060351)
tau2.profile.grid <- c(2.21979442867967, 1.83679414375062, 1.49530794326639, 1.19533471755522, 
0.937228775807105, 0.721860409344911, 0.550900164135041, 0.426702350193188, 
0.35027845740803, 0.318161608378137, 0.323695752454098, 0.362844322995967, 
0.436759149727103, 0.550102939816622, 0.708036309331584, 0.914124795209754, 
1.16984934544558, 1.47523767910086, 1.82951835731162, 2.23187739094619
)
profiled <- structure(c(14.8282746287728, 24.2936658242837, -0.389623832976634, 
13.2675692313222, 15.3711496791249, -0.136852512122726, 11.6042182066143, 
10.6318553540892, 0.0914574944956418, 9.83624770688184, 8.48649822050813, 
0.159046694090144, 7.97129439960355, 6.76599015170296, 0.178141590643201, 
6.03768728773074, 5.07126341622608, 0.190568659559761, 4.10503351404002, 
3.41814193193251, 0.200954669462529, 2.31434699559714, 1.91602623426263, 
0.207888991398807, 0.893917549215097, 0.736038133511206, 0.214498962099614, 
0.102565723520634, 0.0836795411244916, 0.225696533971726, 0.101855068301511, 
0.082637565721863, 0.232551653860775, 0.879553036469375, 0.718139954535075, 
0.224765494406726, 2.279127093215, 1.88199067076217, 0.211019336398731, 
4.07869640896404, 3.3711148911328, 0.209895402761985, 6.06710073471795, 
4.90367131723742, 0.237256810706468, 8.08886704619978, 6.23327479903051, 
0.297691391282458, 10.0516405259884, 7.39458068497647, 0.359325288911952, 
11.9111865464143, -33.901329567223, -1.35134865500763, 13.6524715472944, 
-0.984661176643724, -14.8651465815172, 15.2756241959227, 0.0515842738410368, 
295.129480139556), dim = c(3L, 20L), dimnames = list(c("LR", 
                                                       "LR.bartlett", "bartlett.C"), NULL))
source('utils.R')
ks.ma0 <- function(params,ks.fun,pmf.derivs.fun,mu,tau2,n.GH=50) {
    ## browser()
    ks <- lapply(1:nrow(params), function(i) {
        ## if(i==16)browser()
        supp.x <- params[i,'x.min']:params[i,'x.max']
        ks.df.study <- ks.study0(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params[i,],mu,tau2,n.GH=n.GH) 
    })
    return(ks)
}
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
stijnen.vi <- read.table('stijnen-table-VI.txt')
x <- stijnen.vi$Y0
T <- stijnen.vi$T0/1000
params <- cbind(T=T,x.min=0,x.max=1000)
mu <- mu.profile.grid[20];tau2=tau2.profile.grid[20]
ks <- ks.ma0(params,ks.fun,pmf.derivs.fun,mu=mu,tau2=tau2,n.GH=200)
k.vals.ma <- sapply(ks,function(df)df$k.val) |> rowSums()
ks.profile <- data.frame(inner=ks[[1]]$inner, outer=ks[[1]]$outer,k.val=k.vals.ma)
bartlett.C.profile <- lawley(ks.profile,2) - lawley(ks.profile,1)
print(bartlett.C.profile)
print(summary(ks.profile$k.val))

source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-poisson.rds')
mu <- 2.80764561060351
tau2 <- 2.23187739094619
n.GH <- 50
params <- c(T = 1.988, x.min = 0, x.max = 1e3)
lawley.terms <- sapply(seq(20,3e2,by=20),function(n.GH) {
    cat('.')
    ks.df <- ks.study0(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params,mu,tau2,n.GH)
    ## c(lawley(ks.df,2),lawley(ks.df,1))
    lawley(ks.df,2)
})
lawley.terms

beta <- c(tau2,mu)
## if(is.null(GH)) {
GH <- fastGHQuad::gaussHermiteData(n.GH)
GH$x <- sqrt(2)*GH$x
GH$w <- GH$w / sqrt(pi)
names(GH) <- c('theta.grid','w')
## }
pmf.derivs.unit <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
    ## browser()
    supp.x <-  params['x.min']:params['x.max']
    units.mat <- sapply(GH$theta.grid,function(thetai)pmf.deriv.fun(beta=beta,params=params,theta=thetai,x=supp.x))
    colSums(GH$w * t(units.mat))
})
## print(range(unlist(pmf.derivs.unit)))
## if(max(abs(unlist(pmf.derivs.unit))) > maxval) maxval <<- max(abs(unlist(pmf.derivs.unit)))
k.vals <- sapply(seq_along(ks.fun), function(i) {
    k.fun <- ks.fun[[i]]
    args <- attr(k.fun,'args')
    k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
    k.val <- sum( k.fun.unit )
})
k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
## k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
## a <- pmf.derivs.unit[['pmfSum00']]




n.GHs <- seq(20,2e2)
by.n.GH <- sapply(
beta <- c(tau2,mu)
## if(is.null(GH)) {
GH <- fastGHQuad::gaussHermiteData(n.GH)
GH$x <- sqrt(2)*GH$x
GH$w <- GH$w / sqrt(pi)
names(GH) <- c('theta.grid','w')
## }
pmf.derivs.unit <- lapply(pmf.derivs.fun, function(pmf.deriv.fun) {
    ## browser()
    supp.x <-  params['x.min']:params['x.max']
    units.mat <- sapply(GH$theta.grid,function(thetai)pmf.deriv.fun(beta=beta,params=params,theta=thetai,x=supp.x))
    colSums(GH$w * t(units.mat))
})
## print(range(unlist(pmf.derivs.unit)))
## if(max(abs(unlist(pmf.derivs.unit))) > maxval) maxval <<- max(abs(unlist(pmf.derivs.unit)))
k.vals <- sapply(seq_along(ks.fun), function(i) {
    k.fun <- ks.fun[[i]]
    args <- attr(k.fun,'args')
    k.fun.unit <- do.call(k.fun, pmf.derivs.unit[args])
    k.val <- sum( k.fun.unit )
})
k.orders <- sapply(ks.fun,function(k.fun)attr(k.fun,'order'))
## k.df <- data.frame(inner=sapply(k.orders,`[`,1),outer=sapply(k.orders,`[`,2),k.val=k.vals)
ks.df <- data.frame(inner=k.orders[1,],outer=k.orders[2,],k.val=k.vals)
## a <- pmf.derivs.unit[['pmfSum00']]




REFACTORED

## confint.LR <- function(fit,logtau2.bounds=c(-5,5),grid=list(offset=1,size=20),alpha=.05,n.GH=fit$n.GH) {
##     ## browser()
##     GH <- fastGHQuad::gaussHermiteData(n.GH)
##     GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
##     ## names(GH) <- c('theta.grid','w')
##     x <- fit$x
##     params <- fit$params
##     cond.pmf <- fit$cond.pmf
##     mu.mle.full <- fit$mle.full['mu']; tau2.mle.full <- fit$mle.full['tau2'];
##     loglik <- fit$loglik.fun
##     loglik.profile <- function(mu) {
##         logtau2.mle.profile <- optimize(function(logtau2)loglik(mu,exp(logtau2),x,params,GH),interval=logtau2.bounds,maximum=TRUE)$max
##         tau2.mle.profile <- exp(logtau2.mle.profile)
##         return(tau2.mle.profile)
##     }
##     mu.profile.grid <- seq(mu.mle.full-grid$offset*sqrt(tau2.mle.full),mu.mle.full+grid$offset*sqrt(tau2.mle.full),len=grid$size)
##     tau2.profile.grid <- sapply(mu.profile.grid,function(mu)loglik.profile(mu))
##     ## add support for parallelization
##     profiled <- sapply(seq_along(mu.profile.grid),function(i) {
##         cat('.')
##         LR.profile <- 2 * (loglik(mu.mle.full,tau2.mle.full,x,params,GH) - loglik(mu.profile.grid[i],tau2.profile.grid[i],x,params,GH))
##         ks.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.profile.grid[i],tau2=tau2.profile.grid[i],n.GH=n.GH)
##         bartlett.C.profile <- lawley(ks.profile,2) - lawley(ks.profile,1)
##         inf.stats <- with(list(k=function(inner,outer=0){k.df <- ks.profile; k.df$k.val[k.df$inner==paste(sort(inner),collapse='') & k.df$outer==paste(sort(outer),collapse='')]}),{
##             information <- -matrix(apply(expand.grid(1:2,1:2),1,k),nrow=2)
##             c(det=det(information),min.eigen=min(abs(eigen(information)$val)))
##         })
##         return(c(LR=LR.profile,LR.bartlett=LR.profile / (1+bartlett.C.profile),bartlett.C=bartlett.C.profile,inf.stats))
##     })
##     ## browser()
##     plot(mu.profile.grid,profiled['LR',],ylim=range(profiled[c('LR','LR.bartlett'),]),type='o')
##     lines(mu.profile.grid,profiled['LR.bartlett',],lty=2,type='o')
##     q <- qchisq(1-alpha,df=1)
##     abline(h=q)
##     CIs <- apply( profiled[c('LR','LR.bartlett'),] - q, 1, function(vals) {
##         left.idx <- which(diff(sign(vals))!=0)
##         stopifnot(length(left.idx)==2)
##         stopifnot(left.idx[2]<length(vals))
##         right.idx <- left.idx+1 
##         left.wt <- abs(vals[right.idx]/(vals[right.idx]-vals[left.idx]))
##         right.wt <- abs(vals[left.idx]/(vals[right.idx]-vals[left.idx]))
##         left.wt*mu.profile.grid[left.idx] + right.wt*mu.profile.grid[right.idx]
##     } )
##     return(CIs)
## }




DEBUGGING BINOMIAL-PROP FIT TO NEIL-WEISE. N.GH WASNT HIGH ENOUGH, DETECTING USING INFORMATION MATRIX WHICH ONLY STABILIZED FOR N.GH \APPROX 750

source('utils.R')
ks.fun <- readRDS('ks-fun.rds')
niel.weise <- read.table('niel-weise.txt')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-prop.rds')
x <- niel.weise$Y0
ms <- niel.weise$n0
params <- cbind(m=ms,x.min=0,x.max=ms)
n.GH <- 1e2
dets <- sapply(seq(50,1e3,by=50), function(n.GH) {
    pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-prop.rds')
    mu <- -2.07341681000221
    tau2 <- 3.13080589105322
    ks.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu,tau2=tau2,n.GH=n.GH)
    bartlett.C.profile <- lawley(ks.profile,2) - lawley(ks.profile,1)
    inf.stats <- with(list(k=function(inner,outer=0){k.df <- ks.profile; k.df$k.val[k.df$inner==paste(sort(inner),collapse='') & k.df$outer==paste(sort(outer),collapse='')]}),{
        information <- -matrix(apply(expand.grid(1:2,1:2),1,k),nrow=2)
        c(det=det(information),min.eigen=min(abs(eigen(information)$val)))
    })
    inf.stats['det']
})
print(dets)




SUPERSEDED



diagnostics <- sapply(profile.fit$ks.profile,function(k.df) {
    k=function(inner,outer=0)k.df$k.val[k.df$inner==paste(sort(inner),collapse='') & k.df$outer==paste(sort(outer),collapse='')]
    information <- -matrix(apply(expand.grid(1:2,1:2),1,k),nrow=2)
    c(det=det(information),min.eigen=min(abs(eigen(information)$val)))
})
diagnostics
lines(profile.fit$mu.grid,diagnostics['det',])
lines(profile.fit$mu.grid,profile.fit$profile['bartlett.C',])






n.GH <- 1e2
dets <- sapply(seq(50,1e3,by=50), function(n.GH) {
    pmf.derivs.fun <- readRDS('pmf-derivs-fun-binomial-prop.rds')
    mu <- -2.07341681000221
    tau2 <- 3.13080589105322
    ks.profile <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu,tau2=tau2,n.GH=n.GH)
    bartlett.C.profile <- lawley(ks.profile,2) - lawley(ks.profile,1)
    inf.stats <- with(list(k=function(inner,outer=0){k.df <- ks.profile; k.df$k.val[k.df$inner==paste(sort(inner),collapse='') & k.df$outer==paste(sort(outer),collapse='')]}),{
        information <- -matrix(apply(expand.grid(1:2,1:2),1,k),nrow=2)
        c(det=det(information),min.eigen=min(abs(eigen(information)$val)))
    })
    inf.stats['det']
})
print(dets)
