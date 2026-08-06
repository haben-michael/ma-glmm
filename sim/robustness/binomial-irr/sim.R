args <- commandArgs(trailingOnly=TRUE)
skew <- as.numeric(args[1])
kurt <- as.numeric(args[2])

source('../../../utils.R')
ns <- 5:15

sas.params <- sas.invert.params(skew,kurt)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])

## set.seed(1)#!!!!!!!!!
## re.distr <- rnorm #!!!!

mu.null <- -.476
tau2.null <- .35^2
ks.fun <- readRDS('../../../ks-fun.rds')
pmf.derivs.fun <- readRDS('../../../pmf-derivs-fun-binomial-irr.rds')
## stijnen.vi <- read.table('../../../stijnen-table-VI.txt')
## all.Ts <- sample(stijnen.vi$T0,max(ns),replace=TRUE)
## all.Ts <- all.Ts / 1000
cond.pmf <- cond.pmf.binomial.irr
n.GH <- 50
stijnen.6 <- read.table('../../../stijnen-table-VI.txt')
x <- x.T <- stijnen.6$Y1
x.C <- stijnen.6$Y0
T.t <- stijnen.6$T1/1000
T.c <- stijnen.6$T0/1000
s <- x.C+x.T
params <- cbind(T.t=T.t,T.c=T.c,s=s,x.min=0,x.max=s)
all.params <- params[sample(nrow(params),max(ns),replace=TRUE),]

by.n <- sapply(ns, function(n) {
## repeat{
## by.n <- parallel::mclapply(ns,mc.cores=4,FUN= function(n) {
    cat('.')
    ## Ts <- all.Ts[1:n]
    ## bartletts <- parallel::mclapply(seq(50,5e2,by=50),mc.cores=5,FUN=function(n.GH) {
        ## cat('.')
    params <- all.params[1:n,]
    ks.df <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.null,tau2=tau2.null,n.GH=n.GH)$ks.df
    bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
## })
    stats <- replicate(round(n^3), {
        x <- dgp.binomial.irr(n,mu.null,tau2.null,params,re.distr=re.distr)
        LR <- LR.gaussian(x,params,mu.null,cond.pmf=cond.pmf,n.GH=n.GH)$LR
        LR.bartlett <- LR / (1+bartlett.C)
        c(LR=LR,LR.bartlett=LR.bartlett)
    })
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
    means <- rowMeans(stats-1)
    c(ks=ks,means=means)
})
by.n <- simplify2array(by.n)
## print(Sys.time() - tic)
## plot(ns,by.n['LR',],type='l',ylim=range(by.n))
## lines(ns,by.n['LR.bartlett',],col=2)
## apply(by.n[1:2,],1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()

filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(skew=skew,kurt=kurt,ns=ns,by.n=by.n,file=filename)

