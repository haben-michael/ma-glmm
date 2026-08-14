args <- commandArgs(trailingOnly=TRUE)
skew <- as.numeric(args[1])
kurt <- as.numeric(args[2])
tau2.null <- as.numeric(args[3])

source('../../../utils.R')
ns <- 5:15

sas.params <- sas.invert.params(skew,kurt)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])

# set.seed(13)#!!!!!!!!!
## re.distr <- rnorm #!!!!

mu.null <- -.476
## tau2.null <- .35^2
## mu.null <- 0#-.476
## tau2.null <- 1/2#.35^2
ks.fun <- readRDS('../../../ks-fun.rds')
pmf.derivs.fun <- readRDS('../../../pmf-derivs-fun-binomial-irr.rds')
## stijnen.vi <- read.table('../../../stijnen-table-VI.txt')
## all.Ts <- sample(stijnen.vi$T0,max(ns),replace=TRUE)
## all.Ts <- all.Ts / 1000
cond.pmf <- cond.pmf.binomial.irr
n.GH <- 2e1
stijnen.6 <- read.table('../../../stijnen-table-VI.txt')
stijnen.params <- with(stijnen.6, cbind(T.t=T1,T.c=T0,s=Y0+Y1,x.min=0,x.max=Y0+Y1))
mc.idx <- sample(nrow(stijnen.params),max(ns),replace=TRUE)
mc.params <- stijnen.params[mc.idx,]
## mc.mles <- sapply(ns, function(n)
##     LR.gaussian(stijnen.6$Y1[mc.idx[1:n]],params=stijnen.params[mc.idx[1:n],,drop=FALSE],mu.null=0,cond.pmf=cond.pmf,logtau2.bounds=c(-15,15),n.GH=n.GH,control=list(maxit=1e4))$mle.full
##     ) |> t()
## all.params <- cbind(s=rep(30,max(ns)),T.t=1,T.c=1.5,x.min=0,x.max=30)
## all.params <- do.call(rbind,replicate(max(ns),apply(stijnen.params,2,median),simplify=FALSE))

## by.n <- parallel::mclapply(seq_along(ns),mc.cores=6,FUN= function(i) {
by.n <- sapply(seq_along(ns), function(i) {
    ## repeat{
    cat('.')
    ## Ts <- all.Ts[1:n]
    ## bartletts <- parallel::mclapply(seq(50,5e2,by=50),mc.cores=5,FUN=function(n.GH) {
    ## cat('.')
    ## pairs <- sapply(1:11, function(i) {
    n <- ns[i]
    params <- mc.params[1:n,]
    ## mu.null <-mc.mles[i,'mu']; tau2.null <- mc.mles[i,'tau2']
    ## bartlett.Cs <- sapply(c(20,200), function(n.GH) {
    ks.df <- ks.ma(params,ks.fun,pmf.derivs.fun,mu=mu.null,tau2=tau2.null,n.GH=n.GH)$ks.df
    bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
    ## })
   ## })
    ## })
    stats <- replicate(round(n^4), {
        x <- dgp.binomial.irr(n,mu.null,tau2.null,params,re.distr=re.distr)
        ## idx.cc <- x==0L | x==params[,'s']        
        ## x.T.nn <- x; x.C.nn <- params[,'s']-x;
        ## x.T.nn[idx.cc] <- x.T.nn[idx.cc]+.5; x.C.nn[idx.cc] <- x.C.nn[idx.cc]+.5
        ## y.nn <- log(x.T.nn/params[,'T.t']) - log(x.C.nn/params[,'T.c'])
        ## sigma2.nn <- 1/x.T.nn + 1/x.C.nn
        ## w.fe <- 1/sigma2.nn
        ## mu.fe <- y.nn %*% w.fe / sum(w.fe)    
        ## Q <- w.fe%*%(y.nn-c(mu.fe))^2
        ## tau2.DL <- max(0L, (Q - (length(y.nn)-1)) / (sum(w.fe)-sum(w.fe^2)/sum(w.fe)))
        ## if(is.nan(tau2.DL))browser()
        ## tau2.start <- if (tau2.DL>0L) tau2.DL else 1e-3
        ## mu.start <- mu.fe
        LR <- LR.gaussian(x,params,mu.null,cond.pmf=cond.pmf,logtau2.bounds=c(-15,15),n.GH=n.GH,par.start=par.start.binomial.irr(x,params),control=list(maxit=1e4))$LR
        ## LR <- LR.gaussian(x,params,mu.null,cond.pmf=cond.pmf,n.GH=n.GH)$LR
        ## LR <- LR.gaussian(x,params,mu.null,cond.pmf=cond.pmf,logtau2.bounds=c(-15,15),n.GH=n.GH,control=list(maxit=1e4))$LR
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
## plot(ns,by.n['ks.LR',],type='l',ylim=range(by.n))
## lines(ns,by.n['ks.LR.bartlett',],col=2)
## apply(by.n[1:2,],1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()

filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(tau2.null=tau2.null,skew=skew,kurt=kurt,ns=ns,by.n=by.n,file=filename)

# % while true; do
#   {
#     parallel --jobs -2 Rscript sim.R ::: $(seq -.25 0.05 .25) ::: 0 ::: $(seq .1 .2 .5)
#   }
# done
