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
mu.null <- -1.353
tau2.null <- .83^2
model.str <- 'hn'
n.GH <- 2e1
stijnen.1 <- read.table('../../stijnen-table-I.txt')
stijnen.1 <- stijnen.1[-15,] # get rid of uninformative study
stijnen.params <- with(stijnen.1, cbind(m.T=n1,m.C=n0,s=Y0+Y1))
stijnen.params <- with(as.data.frame(stijnen.params), cbind(stijnen.params,y.min=pmax(0,s-m.C),y.max=pmin(s,m.T)) )                  
mc.idx <- sample(nrow(stijnen.params),max(ns),replace=TRUE)
mc.params <- stijnen.params[mc.idx,]
model0 <- model(model.str)
quad <- make.quad(model0$moment.map,n.GH=n.GH,adaptive=TRUE,adaptive.tol=1e-10)
start <- Sys.time()
by.n <- parallel::mclapply(seq_along(ns),mc.cores=6,FUN= function(i) {
## by.n <- sapply(seq_along(ns), function(i) {
    ## repeat{
    cat('.')
    n <- ns[i]
    params <- mc.params[1:n,]
    bartlett.C <- make.bartlett.C(params,quad,model0$pmf.derivs.fun)
    bartlett.C0 <- bartlett.C(mu.null,tau2.null)$bartlett.C
    stats <- replicate(round(n^4), {        
        y <- model0$dgp(n,mu.null,tau2.null,params,re.distr=re.distr)
        fit <- meta.glmm(y,params,model.str,n.GH,adaptive=FALSE,logtau2.bounds=c(-20,3))
        LR <- fit$profile.LR(mu.null)$LR
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
