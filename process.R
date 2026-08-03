
## process  rate sim (IID data, using R with gauss-hermite)
save.path <- './sim/260714/'
save.path <- './sim/260722/' #poisson
## save.path <- './sim/260723/'
save.path <- './sim/260724/'
save.path <- './sim/260728/'
## png.file <- './figs/gaussian.png'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    load(paste0(save.path,file))
    df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
    do.call(rbind,df.lst)
    ## apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
    ## apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
})
long <- do.call(rbind,by.sim)
long <- aggregate(val ~ ., median,data=long)
plot(val ~ n, ylim=c(0,max(val)),data=subset(long,subset=stat=='ks.LR',select=c(n,val)),type='l')
lines(val ~ n, data=subset(long,subset=stat=='ks.LR.bartlett',select=c(n,val)),col=2)
## plot(val ~ n, ylim=c(0,max(val)),data=subset(long,subset=stat=='means.LR',select=c(n,val)),type='l')
## lines(abs(val) ~ n, data=subset(long,subset=stat=='means.LR.bartlett',select=c(n,val)),col=2)
stats <- split(long,long$stat)
sapply(stats,function(df)coef(lm(log(abs(df$val))~log(df$n)))[2])


## robustness sims

## compute rates, then take median rate
save.path <- './sim/robustness/'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    load(paste0(save.path,file))
    df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
    df <- do.call(rbind,df.lst)
    df <- cbind(df,skew=skew,kurt=kurt)
    ## apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
    rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
    data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
})
long <- do.call(rbind,by.sim)
long <- aggregate(rate ~ ., median,data=long)
long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
df <- long
subset(long,stat=='ks.LR'|stat=='ks.LR.bartlett')
subset(long,stat=='ks.LR')

## compute median stats, then take rates. looks much smoother this way.
save.path <- './sim/robustness/'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    load(paste0(save.path,file))
    df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
    df <- do.call(rbind,df.lst)
    df <- cbind(df,skew=skew,kurt=kurt)
    ## apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
    ## rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
    ## data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
})

long <- do.call(rbind,by.sim)
long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
long <- aggregate(val ~ ., median,data=long)
by.stats <- split(long,long$stat)
## sapply(stats,function(df)coef(lm(log(abs(df$val))~log(df$n)))[2])
by.stats <- lapply(by.stats, function(df)
    rates <- sapply(split(df,list(df$skew,df$kurt)), function(smalldf) {
    if(nrow(smalldf)==0)return(c(NA,NA,NA))
    ## rate <- if(nrow(smalldf)!=0) {
    rate <- with(smalldf, coef(lm(log(abs(val))~log(n)))[[2]])
    ## } else {
        ## NA
    ## }
    c(skew=smalldf$skew[1],kurt=smalldf$kurt[1],rate=rate)
}) |> t()
)
df <- rbind(data.frame(stat='ks.LR',by.stats[['ks.LR']]), data.frame(stat='ks.LR.bartlett',by.stats[['ks.LR.bartlett']]))
rownames(df) <- NULL
df <- na.omit(df)

library(lattice)
levelplot(rate ~ skew * kurt | stat, data = df,
          contour = TRUE,      
          labels = TRUE,       
          pretty = TRUE,       
          region = TRUE,       
          col.regions = gray.colors(100, start = 0.9, end = 0.1),
          layout = c(2, 1),    
          xlab = "Skewness", 
          ylab = "Excess Kurtosis",
          panel = function(...) {
              panel.levelplot(...)
              panel.points(0, 0, pch = 4, col = "red", cex = 1.5) 
          })

subdf <- subset(df,stat=='ks.LR.bartlett' & abs(skew)<=.4 & abs(kurt)<=.4)
matrix(subdf$rate,nrow=length(unique(subdf$skew)))






## skew and kurt separately instead of cartesian product

## process  rate sim (IID data, using R with gauss-hermite)
## save.path <- './sim/robustness/'
save.path <- './sim/robustness/poisson/'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    load(paste0(save.path,file))
    df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
    df <- do.call(rbind,df.lst)
    df <- cbind(df,skew=skew,kurt=kurt)
    ## apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
    ## rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
    ## data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
})
long <- do.call(rbind,by.sim)
long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
long <- aggregate(val ~ ., median,data=long)
by.stats <- split(long,long$stat)
## sapply(stats,function(df)coef(lm(log(abs(df$val))~log(df$n)))[2])
by.stats <- lapply(by.stats, function(df) 
    ## rates <- sapply(split(df,list(df$skew,df$kurt)), function(smalldf)
    rates <- sapply(split(df,paste0(df$skew,df$kurt)), function(smalldf) {
        if(nrow(smalldf)==0)return(c(NA,NA,NA))
        ## rate <- if(nrow(smalldf)!=0) {
        rate <- with(smalldf, coef(lm(log(abs(val))~log(n)))[[2]])
        ## } else {
        ## NA
        ## }
        c(skew=smalldf$skew[1],kurt=smalldf$kurt[1],rate=rate)
    }) |> t()    
    )
df <- rbind(data.frame(stat='ks.LR',by.stats[['ks.LR']]), data.frame(stat='ks.LR.bartlett',by.stats[['ks.LR.bartlett']]))
rownames(df) <- NULL
## df <- na.omit(df)

df <- df[order(df$skew,df$kurt),]
with(subset(df,kurt==0 & stat=='ks.LR'),  plot(rate ~ skew,type='l',ylim=range(df$rate)))
with(subset(df,kurt==0 & stat=='ks.LR.bartlett'),  lines(rate ~ skew,col=2))

df <- df[order(df$skew,df$kurt),]
with(subset(df,skew==0 & stat=='ks.LR'),  plot(rate ~ kurt,type='l',ylim=range(df$rate)))
with(subset(df,skew==0 & stat=='ks.LR.bartlett'),  lines(rate ~ kurt,col=2))
