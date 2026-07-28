source('../../utils.R')
## check bartlett rates        
## set.seed(1)
mu.true <- mu.null <- 0
tau2.true <- tau2.null <- 1/2
## n <- 10
ns <- 5:15
## bartlett.C <- 1.36433398617429 # hyperg
## bartlett.C <- 1.17164232232645 # poisson
bartlett.C <- 1.36297164231076 # binomial-IRR
## tic <- Sys.time()
by.n <- sapply(ns, function(n) {
## repeat{
## by.n <- parallel::mclapply(ns,mc.cores=4,FUN= function(n) {
    cat('.')
    s <- rep(30,n)
    T.t <- rep(1,n)
    T.c <- rep(1.5,n)
    params <- cbind(s=s,T.t=T.t,T.c=T.c)
    stats <- replicate(round(n^4), {
        x <- dgp.binomial.irr(n,mu.null,tau2.null,params)
        LR <- LR.gaussian(x,params,mu.null,cond.pmf.binomial.irr)$LR
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

