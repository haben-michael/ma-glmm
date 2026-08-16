args <- commandArgs(trailingOnly=TRUE)
skew <- as.numeric(args[1])
kurt <- as.numeric(args[2])
tau2.null <- as.numeric(args[3])

source('../../../utils.R')
ns <- 5:15

sas.params <- sas.invert.params(skew,kurt)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])

RNGkind("L'Ecuyer-CMRG")#!!!!!!
set.seed(1)#!!!!!!!!!
re.distr <- rnorm #!!!!
tau2.null <- 1/2
## mu.null <- 0

mu.null <- -1.353
## tau2.null <- .83^2
# ks.fun <- readRDS('../../../ks-fun.rds')
# pmf.derivs.fun <- readRDS('../../../pmf-derivs-fun-hyperg.rds')
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun <- readRDS('pmf-derivs-fun-hyperg.rds')
cond.pmf <- cond.pmf.hyperg
dgp <- dgp.hyperg
n.GH <- 4e1
stijnen.1 <- read.table('stijnen-table-I.txt')
stijnen.1 <- stijnen.1[-15,] # get rid of uninformative study
stijnen.params <- with(stijnen.1, cbind(m.T=n1,m.C=n0,s=Y0+Y1))
stijnen.params <- with(as.data.frame(stijnen.params), cbind(stijnen.params,x.min=pmax(0,s-m.C),x.max=pmin(s,m.T)) )                  
mc.idx <- sample(nrow(stijnen.params),max(ns),replace=TRUE)
mc.params <- stijnen.params[mc.idx,]
## mc.mles <- sapply(ns, function(n)
    ## LR.gaussian(stijnen.1$Y1[mc.idx[1:n]],params=stijnen.params[mc.idx[1:n],,drop=FALSE],mu.null=0,cond.pmf=cond.pmf,logtau2.bounds=c(-15,15),n.GH=n.GH,control=list(maxit=1e4))$mle.full
    ## ) |> t()
## all.params <- cbind(s=rep(30,max(ns)),T.t=1,T.c=1.5,x.min=0,x.max=30)
## all.params <- do.call(rbind,replicate(max(ns),apply(stijnen.params,2,median),simplify=FALSE))

start <- Sys.time()
# by.n <- parallel::mclapply(seq_along(ns),mc.cores=parallel::detectCores()-1,FUN= function(i) {
 by.n <- sapply(seq_along(ns), function(i) {
    cat('.')
    ## Ts <- all.Ts[1:n]
    ## bartletts <- parallel::mclapply(seq(50,5e2,by=50),mc.cores=5,FUN=function(n.GH) {
    ## cat('.')
    ## pairs <- sapply(1:11, function(i) {
    n <- ns[i]
    params <- mc.params[1:n,]
    paramsi <- cbind(m.T = 20, m.C = 15, s = 10, x.min = 0, x.max = 10)#!!!!!!
    paramsi <- cbind(m.T = 100, m.C = 100, s = 10, x.min = 0, x.max = 10)#!!!!!!
    params <- paramsi[rep(1,n),]
    ## params <- mc.params[rep(1,n),] #!!!!!!!!!!!!!!!
    ## mu.null <-mc.mles[i,'mu']; tau2.null <- mc.mles[i,'tau2']
    ## bartlett.Cs <- sapply(c(20,200), function(n.GH) {
    ks.df <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.null,tau2=tau2.null,n.GH=n.GH)$ks.df
    bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
    stats <- replicate(round(n^3), {
        x <- dgp(n,mu.null,tau2.null,params,re.distr=re.distr)
        ## LR <- LR.gaussian(x,params,mu.null,cond.pmf=cond.pmf,n.GH=n.GH)$LR
        LR.fit <- LR.gaussian(x,params,mu.null,cond.pmf=cond.pmf,logtau2.bounds=c(-20,3),n.GH=n.GH,par.start=par.start.hyperg(x,params),control=list(maxit=1e4))
        LR <- LR.fit$LR
        LR.bartlett <- LR / (1+bartlett.C)
        c(LR=LR,LR.bartlett=LR.bartlett)#,LR.fit$mle.full)
    })
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
    means <- rowMeans(stats-1)
    c(ks=ks,means=means)
})
by.n <- simplify2array(by.n)
print(Sys.time() - start)
plot(ns,by.n['ks.LR',],type='l',ylim=range(by.n[1:2,]))
lines(ns,by.n['ks.LR.bartlett',],col=2)
apply(by.n[1:2,],1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()

filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(tau2.null=tau2.null,skew=skew,kurt=kurt,ns=ns,by.n=by.n,file=filename)

# % while true; do
#   {
#     parallel --jobs -2 Rscript sim.R ::: $(seq -.25 0.05 .25) ::: 0 ::: $(seq .1 .2 .5)
#   }
# done

n <- 30
ks.df <- ks.study(pmf.derivs.fun,ks.fun,params[1,],mu.null,tau2.null,n.GH)$ks.df
k.22 <- ks.df[ks.df[,'inner']=='22'  & ks.df[,'outer']=='0','k.val']
k.11 <- ks.df[ks.df[,'inner']=='11'  & ks.df[,'outer']=='0','k.val']
k.12 <- ks.df[ks.df[,'inner']=='12'  & ks.df[,'outer']=='0','k.val']
inf <- -matrix(c(k.11,k.12,k.12,k.22),nrow=2)
vcov.hat <- solve(inf)
beta <- c(tau2.null,mu.null)
z.stats <- replicate(1e2, {
    x <- dgp.hyperg(n,mu.null,tau2.null,params)
    mles <- LR.gaussian(x,params,mu.null,cond.pmf.hyperg,n.GH=n.GH)$mle.full
    (rev(mles)-beta) / sqrt(diag(vcov.hat)/n)
})
op <- par(mfrow=c(1,2))
qqnorm(z.stats['mu',]); abline(0,1)
qqnorm(z.stats['tau2',]); abline(0,1)
par(op)


cond.pmf <- function(x,theta,mu,tau2,params) {
    m.T <- params[,'m.T']
    m.C <- params[,'m.C']
    s <- params[,'s']
    odds <- exp(mu+sqrt(tau2)*theta)
    out <- .cond.pmf.hyperg(x=x,m1=m.T,m2=m.C,n=s,odds=odds)
    ## print(out)
    return(out)
}
params <- structure(c(100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 
10, 10, 10, 10, 10, 0, 0, 0, 0, 0, 10, 10, 10, 10, 10), dim = c(5L, 
5L), dimnames = list(NULL, c("m.T", "m.C", "s", "x.min", "x.max"
)))
x <- c(0L, 1L, 0L, 9L, 0L)
mu <- c(mu = -3.78862407519942)
tau2 <- c(tau2 = 23.7842558474395)
n.GH <- 5e2
GH <- fastGHQuad::gaussHermiteData(n.GH)
GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
logtau2.bounds=c(-20,5)
loglik.unit <- function(mu,tau2,xi,paramsi,GH) {
    ## browser()
    ## print(c(mu,tau2,xi,paramsi))
    ## if(tau2>1e3)browser()
    if(log(tau2)<logtau2.bounds[1] || log(tau2)>logtau2.bounds[2]) {
        ## warning('tau2 boundary  hit, returning minimum negative')
        return(-1e10)#     return(-.Machine$double.xmax)
    }
    likelihood <- sum(GH$w*cond.pmf(xi,GH$x,mu,tau2,paramsi))
    ## print(c(mu,tau2,log(likelihood),paramsi))
    ## print(tau2)
    return(log(likelihood))
}
loglik <- function(mu,tau2,x,params,GH) {
    out <- sum(sapply(seq_along(x),function(i)loglik.unit(mu,tau2,x[i],params[i,,drop=FALSE],GH)))
    print(c(mu,tau2,out))
    ## if(abs(out)>1e6)browser()
    ## print(out)        
    ## if(is.nan(out))browser()
    if(is.nan(out))out <- -1e10#-.Machine$double.xmax
    return(out)
}
par.start <- c(mu=-1.402138,tau2=  5.232188 )
fit.full <- optim(par=c(par.start['mu'],log(par.start['tau2'])),fn=function(par)-loglik(par[1],exp(par[2]),x,params,GH),method='BFGS',control=list(maxit=1e4,trace=5))



lapply(seq_along(x), function(i) {
    p <- cond.pmf(
        x[i], GH$x,
        par.start["mu"], par.start["tau2"],
        params[i,,drop=FALSE]
    )

    bad <- which(!is.finite(p))

    if (!length(bad))
        return(NULL)

    data.frame(
        i = i,
        x = x[i],
        theta = GH$x[bad],
        weight = GH$w[bad],
        logweight = log(GH$w[bad]),
        eta = par.start["mu"] +
              sqrt(par.start["tau2"]) * GH$x[bad],
        pmf = p[bad]
    )
})




f <- function(x)exp(x)
integrate(function(x)f(x)*dnorm(x),-Inf,Inf)$val
integrate(function(x)f(x)*dnorm(x),-1e2,1e2)$val
n.GH <- 5e2
GH <- fastGHQuad::gaussHermiteData(n.GH)
GH$x <- sqrt(2)*GH$x; GH$w <- GH$w/sqrt(pi)
sum(f(GH$x)*GH$w)


