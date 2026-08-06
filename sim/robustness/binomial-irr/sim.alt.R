args <- commandArgs(trailingOnly=TRUE)
skew <- as.numeric(args[1])
kurt <- as.numeric(args[2])

source('../../../utils.R')
ns <- 5:15

sas.params <- sas.invert.params(skew,kurt)
re.distr <- function(n)rsas(n,eps=sas.params['eps'],delta=sas.params['delta'])

## set.seed(1)#!!!!!!!!!


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
stijnen.6.params <- cbind(T.t=T.t,T.c=T.c,s=s,x.min=0,x.max=s)
params <- stijnen.6.params[sample(nrow(stijnen.6.params),max(ns),replace=TRUE),]


x.mat <- replicate(max(ns)^3,dgp.binomial.irr(nrow(params),mu=mu.null,tau2=tau2.null,params=params,re.distr=rnorm)) 


## start <- Sys.time()
LRs <- sapply(seq_along(ns), function(i) {
    n <- ns[i]
    sample <- x.mat[1:n,]
    apply(sample,2,function(x)LR.gaussian(x,params=params[1:n,],mu.null=mu.null,cond.pmf=cond.pmf,n.GH=n.GH)$LR)
})
## print(Sys.time() - start)


bartlett.Cs <-  sapply(seq_along(ns), function(i) {
    n <- ns[i]
    ks.df <- ks.ma(params[1:n,],ks.fun,pmf.derivs.fun,mu=mu.null,tau2=tau2.null,n.GH=n.GH)$ks.df.ma
    bartlett.C <- lawley(ks.df,2) - lawley(ks.df,1)
})

LR.bartletts <- LRs / (1+bartlett.Cs)

ecdfs <- apply(LRs,2,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
ks.LR <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))

ecdfs <- apply(LR.bartletts,2,function(x)ecdf(pchisq(x,df=1,lower=FALSE)))
ks.LR.bartlett <- sapply(ecdfs,function(ecdf)max(abs(ecdf(knots(ecdf)) - knots(ecdf))))

by.n <- rbind(LR=ks.LR,LR.bartlett=ks.LR.bartlett)

## apply(by.n[1:2,],1,function(ks)coef(lm(log(ks) ~ log(ns)))[2]) |> print()
filename <- paste0('save',as.integer(abs(rnorm(1))*1e8),'.RData')
save(skew=skew,kurt=kurt,ns=ns,by.n=by.n,file=filename)

