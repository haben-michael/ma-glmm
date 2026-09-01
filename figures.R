
## poisson convergence matrix

load('pn diagnostic--adaptive--trt 4.RData')
mat <- lapply(CIs[!is.na(CIs)],`[[`,'CIs') |> sapply(`[[`,'LR.bartlett')  |> t() |> round(4) |> apply(1,paste,collapse=',') |> matrix(nrow=length(unique(pars$n.GH)))
print(mat)

matrix.to.latex.table <- function(mat, pars) {
  ## Helper to pad a numeric string to 4 decimal places
  format.val <- function(val.str) {
    num <- as.numeric(val.str)
    return(sprintf("%.4f", num))
  }
  ## Process cells: split by comma, pad numbers, wrap in brackets
  format.cell <- function(cell.str) {
    parts <- unlist(strsplit(cell.str, ","))
    return(paste0("[", format.val(parts[1]), ", ", format.val(parts[2]), "]"))
  }
  ## Apply formatting across the entire matrix
  formatted.mat <- apply(mat, c(1, 2), format.cell)
  ## Extract unique row and column values from pars
  r.vals <- unique(pars$n.GH)
  c.vals <- unique(pars$y.max)
  n.cols <- length(c.vals)
  ## Print LaTeX output
  ## cat("\\begin{table}[htbp]\n\\centering\n")
  cat(sprintf("\\begin{tabular}{c %s}\n\\toprule\n", paste(rep("c", n.cols), collapse = " ")))
  cat(sprintf("& \\multicolumn{%d}{c}{\\textbf{Poisson support cutoff $y_{\\max}$}} \\\\\n", n.cols))
  cat(sprintf("\\cmidrule(lr){2-%d}\n", n.cols + 1))
  cat(sprintf("\\textbf{$n_{\\text{GH}}$} & %s \\\\\n\\midrule\n", paste(sprintf("\\textbf{%s}", c.vals), collapse = " & ")))
  for (i in 1:nrow(formatted.mat)) {
    row.str <- paste(formatted.mat[i, ], collapse = " & ")
    cat(sprintf("%s & %s \\\\\n", r.vals[i], row.str))
  }
    ## cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
      cat("\\bottomrule\n\\end{tabular}\n")
}
sink('figs/pn-table.tex')
matrix.to.latex.table(mat, pars)
sink()





## pn CIs--adaptive
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun.pn <- readRDS('pmf-derivs-fun-poisson.rds')
source('models.R')
source('./meta.glmm/R/meta.glmm.R')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 1e2
y <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
params <- cbind(T=T,y.min=0,y.max=1e4)
n.GH <- c(seq(10,20,by=5),1e2)
y.max <- c(seq(5e2,1.5e3,by=5e2),5e3)#seq(1e2,3e2,by=50)
pars <- expand.grid(n.GH=n.GH,y.max=y.max)
pars <- expand.grid(n.GH=c(20,1e2),y.max=5e3)
CIs <- apply(pars, 1, function(par) {
    print(par)
    n.GH <- par['n.GH']; y.max <- par['y.max']
    params <- cbind(T=T,y.min=0,y.max=y.max)
    fit <- meta.glmm(y,params,'pn',n.GH=n.GH,adaptive=TRUE,logtau2.bounds=c(-20,3))
    CI <- tryCatch({
        fit <- confint.meta.glmm(fit,mc.cores=parallel::detectCores()-1,mu.grid=list(offset=1.5,size=20),plot=FALSE)
        print(fit)
        fit
    },    error=function(e)NA)
    return(CI)
})

## load('pn diagnostic--adaptive--trt 4.RData')
load('pn diagnostic--adaptive--trt 5.RData')
alpha <- .05
mat <- lapply(CIs[!is.na(CIs)],`[[`,'CIs') |> sapply(`[[`,'LR.bartlett')  |> t() |> round(4) |> apply(1,paste,collapse=',') |> matrix(nrow=length(unique(pars$n.GH)))
print(mat)
CIs <- CIs[c(1,6,11,16)]
mu.grid <- CIs[[1]]$mu.grid
## LR <- CIs[[1]]$LR
LR <- sapply(lapply(CIs,`[[`,'profile'),`[`,'LR',)
LR.bartlett <- sapply(lapply(CIs,`[[`,'profile'),`[`,'LR.bartlett',)
matplot(mu.grid,LR,type='l',lty=2,ylim=range(c(LR,LR.bartlett)),col=1,xlab=expression(mu),ylab='profile likelihood')
abline(h=qchisq(1-alpha,df=1))
matplot(mu.grid,LR.bartlett,type='l',lty=1,add=TRUE,col=1)





## nonadaptive quadrature
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun.pn <- readRDS('pmf-derivs-fun-poisson.rds')
source('models.R')
source('./meta.glmm/R/meta.glmm.R')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 1e2
y <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
params <- cbind(T=T,y.min=0,y.max=1e4)
n.GH <- seq(100,500,by=200)
y.max <- seq(1e3,1.3e4,by=3e3)
pars <- expand.grid(n.GH=n.GH,y.max=y.max)
CIs <- apply(pars, 1, function(par) {
    print(par)
    n.GH <- par['n.GH']; y.max <- par['y.max']
    params <- cbind(T=T,y.min=0,y.max=y.max)
    fit <- meta.glmm(y,params,'pn',n.GH=n.GH,adaptive=FALSE,logtau2.bounds=c(-20,3))
    CI <- tryCatch({
        fit <- confint.meta.glmm(fit,mc.cores=parallel::detectCores()-1,mu.grid=list(offset=2.5,size=10),plot=FALSE)
        print(fit)
        fit
    },    error=function(e)NA)
    return(CI)
})
## save.image('pn diagnostic--nonadaptive--trt.RData')

load('pn diagnostic--nonadaptive--trt.RData')
sapply(lapply(lapply(CIs,`[[`,'CIs'),`[[`,'LR.bartlett'),function(mat)if(nrow(mat)==2)mat[2,] else mat[1,]) |> round(4) |> apply(2,paste,collapse=',') |> matrix(nrow=length(unique(pars$n.GH)))




alpha <- .05
mu.grid <- CIs[[1]]$mu.grid
## LR <- CIs[[1]]$LR
LR <- sapply(lapply(CIs[pars$y.max==13000],`[[`,'profile'),`[`,'LR',)
LR.bartlett <- sapply(lapply(CIs[pars$y.max==13000],`[[`,'profile'),`[`,'LR.bartlett',)
matplot(mu.grid,LR,type='l',lty=2,ylim=range(c(LR,LR.bartlett)),col=1,xlab=expression(mu),ylab='likelihood')
abline(h=qchisq(1-alpha,df=1))
matplot(mu.grid,LR.bartlett,type='l',lty=1,add=TRUE,col=1)





## diagnostic plot for pmf derivatives and k summands

## pn CIs--adaptive
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun.pn <- readRDS('pmf-derivs-fun-poisson.rds')
source('models.R')
source('./meta.glmm/R/meta.glmm.R')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 1e2
y <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
## params <- cbind(T=T,y.min=0,y.max=1e4)
## n.GH <- c(seq(10,20,by=5),1e2)
## y.max <- c(seq(5e2,1.5e3,by=5e2),5e3)#seq(1e2,3e2,by=50)
n.GH <- 20
y.max <- 3.5e2
params <- cbind(T=T,y.min=0,y.max=y.max)
fit <- meta.glmm(y,params,'pn',n.GH=n.GH,adaptive=TRUE,logtau2.bounds=c(-20,3),debug=TRUE)
## CI <- confint.meta.glmm(fit,mc.cores=parallel::detectCores()-1,mu.grid=list(offset=1.5,size=10),plot=TRUE)
## mle.null <- fit$profile.LR( -1.1726610)$mle.null
bartlett.fit <- fit$bartlett.C(fit$mle['mu'],fit$mle['tau2'])
## diagnostic.k(bartlett.fit,3)

png('./figs/diagnostic-pmf.png', width = 7, height = 5, units = "in", res = 300, pointsize = 12)
diagnostic.pmf(bartlett.fit,3)
dev.off()

png('./figs/diagnostic-k-adaptive.png', width = 9, height = 7, units = "in", res = 300, pointsize = 12)
diagnostic.k(bartlett.fit,3)
dev.off()

## pn CIs--non-adaptive
ks.fun <- readRDS('ks-fun.rds')
pmf.derivs.fun.pn <- readRDS('pmf-derivs-fun-poisson.rds')
source('models.R')
source('./meta.glmm/R/meta.glmm.R')
stijnen.vi <- read.table('stijnen-table-VI.txt')
n.GH <- 1e2
y <- stijnen.vi$Y1
T <- stijnen.vi$T1/1000
## params <- cbind(T=T,y.min=0,y.max=1e4)
## n.GH <- seq(100,500,by=200)
## y.max <- seq(1e3,1.3e4,by=3e3)
## pars <- expand.grid(n.GH=n.GH,y.max=y.max)
## CIs <- apply(pars, 1, function(par) {
## par <- unlist(pars[13,])
## print(par)
## n.GH <- par['n.GH']; y.max <- par['y.max']
n.GH <- 20;y.max <- 3e3
params <- cbind(T=T,y.min=0,y.max=y.max)
fit <- meta.glmm(y,params,'pn',n.GH=n.GH,adaptive=FALSE,logtau2.bounds=c(-20,3),debug=TRUE)
## CI <- confint.meta.glmm(fit,mc.cores=parallel::detectCores()-1,mu.grid=list(offset=2.5,size=10),plot=TRUE)
## mle.null <- fit$profile.LR( -6)$mle.null
bartlett.fit <- fit$bartlett.C(fit$mle['mu'],fit$mle['tau2'])

study.no <- 3
## png('./figs/diagnostic-k-nonadaptive.png')
png("./figs/diagnostic-k-nonadaptive.png", width = 9, height = 7, units = "in", res = 300, pointsize = 12)
vals <- bartlett.fit$debug[[study.no]]$ks.unit[,]
matplot(asinh(vals),type='l',pch=1,col=1,xlab='y',ylab='k',lty=1,yaxt='n')
## ticks <- round(quantile(c(vals)))
## ticks <- range(vals); ticks <- seq(ticks[1],ticks[2],len=4); ticks <- round(ticks,-1)
ticks <- quantile(vals,c(.001,.02,.5,.98,.999)) |> round()
## print(ticks)
axis(2, at = asinh(ticks), labels = ticks)
dev.off()
