args <- commandArgs(trailingOnly=TRUE)
skew <- as.numeric(args[1])
kurt <- as.numeric(args[2])

source('../../utils.R')
ns <- 5:15

sas.params <- sas.invert.params(skew,kurt)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])


mu.null <- -3.5
tau2.null <- .81
ks.fun <- readRDS('../../ks-fun.rds')
pmf.derivs.fun <- readRDS('../../pmf-derivs-fun-binomial-prop.rds')
niel.weise <- read.table('../../niel-weise.txt')
all.ms <- sample(niel.weise$n0,max(ns),replace=TRUE)

by.n <- sapply(ns, function(n) {
## repeat{
## by.n <- parallel::mclapply(ns,mc.cores=4,FUN= function(n) {
    cat('.')
    ms <- all.ms[1:n]
    params <- cbind(m=ms)
    ks <- lapply(ms, function(m) {
        supp.x <- 0:m
        params <- c(m=m)
        ks.df.study <- ks.study(pmf.derivs.fun=pmf.derivs.fun,ks.fun=ks.fun,params=params,supp.x=supp.x,mu.null=mu.null,tau2.null=tau2.null,n.GH=3e2)    
    })
    k.vals.ma <- sapply(ks,function(df)df$k.val) |> rowSums()
    ks.df.ma <- data.frame(inner=ks[[1]]$inner, outer=ks[[1]]$outer,k.val=k.vals.ma)
    bartlett.C <- lawley(ks.df.ma,2) - lawley(ks.df.ma,1)
    stats <- replicate(round(n^4), {
        x <- dgp.binomial.prop(n,mu.null,tau2.null,params,re.distr=re.distr)
        LR <- LR.gaussian(x,params,mu.null,cond.pmf.binomial.prop,n.GH=1e2)$LR
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
## apply(by.n,1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()

filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(skew=skew,kurt=kurt,ns=ns,by.n=by.n,file=filename)

