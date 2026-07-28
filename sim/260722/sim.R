source('../../utils.R')
## check bartlett rates        
## set.seed(1)
mu.null <- -4
tau2.null <- .2
## n <- 10
ns <- 5:15
## bartlett.C <- 1.36433398617429 # hyperg
bartlett.C <- 1.47249385400558 # poisson
## tic <- Sys.time()
by.n <- sapply(ns, function(n) {
## repeat{
## by.n <- parallel::mclapply(ns,mc.cores=4,FUN= function(n) {
    cat('.')
    ## m.T <- rep(20,n)
    ## m.C <- rep(15,n)
    ## s <- rep(10,n)
    ## s <- rep(10,n)
    T <- rep(1e3,n)
    stats <- replicate(round(n^4), {
        theta <- rnorm(n)
        ## T <- 4 + rexp(n)
        lambda <- T*exp(mu.null+sqrt(tau2.null)*theta)
        x <- rpois(n,lambda)
        params <- cbind(T=T)
        LR <- LR.gaussian(x,params,mu.null,cond.pmf.poisson,n.GH=1e2)$LR
        ## LR <- LR.gaussian(x,params,mu.null,loglik.unit.hyperg)$LR
        LR.bartlett <- LR / (1+bartlett.C/n)
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
save(ns=ns,by.n=by.n,file=filename)

