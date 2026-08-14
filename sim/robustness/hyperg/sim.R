args <- commandArgs(trailingOnly=TRUE)
skew <- as.numeric(args[1])
kurt <- as.numeric(args[2])
tau2.null <- as.numeric(args[3])

source('../../../utils.R')
ns <- 5:15

sas.params <- sas.invert.params(skew,kurt)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])

## RNGkind("L'Ecuyer-CMRG")#!!!!!!
## set.seed(1)#!!!!!!!!!
## re.distr <- rnorm #!!!!
## tau2.null <- 3

mu.null <- -1.353
## tau2.null <- .83^2
ks.fun <- readRDS('../../../ks-fun.rds')
pmf.derivs.fun <- readRDS('../../../pmf-derivs-fun-hyperg.rds')
cond.pmf <- cond.pmf.hyperg
dgp <- dgp.hyperg
n.GH <- 4e1
stijnen.1 <- read.table('../../../stijnen-table-I.txt')
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
## by.n <- parallel::mclapply(seq_along(ns),mc.cores=6,FUN= function(i) {
by.n <- sapply(seq_along(ns), function(i) {
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
## plot(ns,by.n['ks.LR',],type='l',ylim=range(by.n[1:2,]))
## lines(ns,by.n['ks.LR.bartlett',],col=2)
## apply(by.n[1:2,],1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()

filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(tau2.null=tau2.null,skew=skew,kurt=kurt,ns=ns,by.n=by.n,file=filename)

# % while true; do
#   {
#     parallel --jobs -2 Rscript sim.R ::: $(seq -.25 0.05 .25) ::: 0 ::: $(seq .1 .2 .5)
#   }
# done

