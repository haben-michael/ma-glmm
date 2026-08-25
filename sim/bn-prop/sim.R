args <- commandArgs(trailingOnly=TRUE)
skew <- as.numeric(args[1])
kurt <- as.numeric(args[2])
tau2.null <- as.numeric(args[3])

ks.fun <- readRDS('../../ks-fun.rds')
pmf.derivs.fun.pn <- readRDS('../../pmf-derivs-fun-poisson.rds')
pmf.derivs.fun.bn.prop <- readRDS('../../pmf-derivs-fun-binomial-prop.rds')
pmf.derivs.fun.bn.irr <- readRDS('../../pmf-derivs-fun-binomial-irr.rds')
pmf.derivs.fun.hn <- readRDS('../../pmf-derivs-fun-hyperg.rds')
source('../../models.R')
source('../../meta.glmm/R/meta.glmm.R')


sas.params <- sas.invert.params(skew,kurt)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])
## RNGkind("L'Ecuyer-CMRG")#!!!!!!
## set.seed(1)#!!!!!!!!!
## re.distr <- rnorm #!!!!
## tau2.null <- 1/2


ns <- 5:15
mu.null <- -3.5
tau2.null <- .81
model.str <- 'bn-prop'
n.GH <- 1e2
stijnen.1 <- read.table('../../stijnen-table-I.txt')
y <- stijnen.1$Y0
ms <- stijnen.1$n0
stijnen.params <- cbind(m=ms,y.min=0,y.max=ms)
mc.idx <- sample(nrow(stijnen.params),max(ns),replace=TRUE)
mc.params <- stijnen.params[mc.idx,]
model0 <- model(model.str)
quad <- make.quad(model0$moment.map,model0$cond.pmf,n.GH=n.GH,adaptive=FALSE,adaptive.tol=1e-5)
start <- Sys.time()
## by.n <- parallel::mclapply(seq_along(ns),mc.cores=6,FUN= function(i) {
by.n <- sapply(seq_along(ns), function(i) {
    ## repeat{
    cat('.')
    n <- ns[i]
    params <- mc.params[1:n,]
    bartlett.C <- make.bartlett.C(params,quad,model0$pmf.derivs.fun)
    bartlett.C0 <- bartlett.C(mu.null,tau2.null)$bartlett.C
    stats <- replicate(round(n^4), {        
        y <- model0$dgp(n,mu.null,tau2.null,params,re.distr=re.distr)
        # fit <- meta.glmm(y,params,model.str,n.GH,adaptive=FALSE,logtau2.bounds=c(-20,3))        
        # LR <- fit$profile.LR(mu.null)$LR
        loglik <- make.loglik(model0$cond.pmf,params,y,quad)
        logtau2.mle.null <- optimize(function(logtau2)loglik(mu.null,exp(logtau2)),interval=c(-20,3),maximum=TRUE)$max
        tau2.mle.null <- exp(logtau2.mle.null)
        par.start0 <- c(mu=mu.null+rnorm(1)/sqrt(n),tau2=tau2.mle.null)
        mle0 <- mle(loglik,par.start0,method='BFGS')$mle
        profile.LR <- make.profile.LR(loglik,mle0,par.start0,logtau2.bounds=c(-20,3))
        LR <- profile.LR(mu.null)$LR
        LR.bartlett <- LR / (1+bartlett.C0)
        c(LR=LR,LR.bartlett=LR.bartlett)
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
