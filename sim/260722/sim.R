source('../../utils.R')
## check bartlett rates        
## set.seed(1)
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
## n <- 10
ns <- 5:15
## bartlett.C <- 1.36433398617429 # hyperg
bartlett.C <- 1.17164232232645 # poisson
## tic <- Sys.time()
## by.n <- sapply(ns, function(n) {
## repeat{
by.n <- parallel::mclapply(ns,mc.cores=4,FUN= function(n) {
    cat('.')
    ## m.T <- rep(20,n)
    ## m.C <- rep(15,n)
    ## s <- rep(10,n)
    ## s <- rep(10,n)
    T <- rep(1,n)
    stats <- replicate(round(n^4), {
        theta <- rnorm(n)
        ## T <- 4 + rexp(n)
        lambda <- T*exp(mu.null+sqrt(tau2.null)*theta)
        x <- rpois(n,lambda)
        params <- cbind(T=T)
        ## theta <- rnorm(n)
        ## odds <- exp(mu.true+sqrt(tau2.true)*theta)
        ## x <- x.T <- sapply(1:n, function(i) BiasedUrn::rFNCHypergeo(1,m1=m.T[i],m2=m.C[i],n=s[i],odds=odds[i]))
        ## params <- cbind(m.T=m.T,m.C=m.C,s=s)
        LR <- LR.gaussian(x,params,mu.null,cond.pmf.poisson)$LR
        ## LR <- LR.gaussian(x,params,mu.null,loglik.unit.hyperg)$LR
        LR.bartlett <- LR / (1+bartlett.C/n)
        c(LR=LR,LR.bartlett=LR.bartlett)
    })
    ecdfs <- apply(stats,1,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
    ks <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))
})
by.n <- simplify2array(by.n)
## print(Sys.time() - tic)
## plot(ns,by.n['LR',],type='l',ylim=range(by.n))
## lines(ns,by.n['LR.bartlett',],col=2)
## apply(by.n,1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()

filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(ns=ns,by.n=by.n,file=filename)

