
## process  rate sim (IID data, using R with gauss-hermite)
save.path <- './sim/260714/'
save.path <- './sim/260722/' #poisson
save.path <- './sim/260723/'
## save.path <- './sim/260724/'
## save.path <- './sim/260728/'
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


## ## robustness sims

## ## compute rates, then take median rate
## save.path <- './sim/robustness/'
## filelist <- dir(save.path)
## filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
## print(length(filelist))
## by.sim <- lapply(filelist, function(file) {
##     load(paste0(save.path,file))
##     df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
##     df <- do.call(rbind,df.lst)
##     df <- cbind(df,skew=skew,kurt=kurt)
##     ## apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
##     rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
##     data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
## })
## long <- do.call(rbind,by.sim)
## long <- aggregate(rate ~ ., median,data=long)
## long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
## df <- long
## subset(long,stat=='ks.LR'|stat=='ks.LR.bartlett')
## subset(long,stat=='ks.LR')

## ## compute median stats, then take rates. looks much smoother this
## ## way. but doesnt make sense a priori--k-s rates across different
## ## designs are being averaged/median'd.
## save.path <- './sim/robustness/'
## filelist <- dir(save.path)
## filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
## print(length(filelist))
## by.sim <- lapply(filelist, function(file) {
##     load(paste0(save.path,file))
##     df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
##     df <- do.call(rbind,df.lst)
##     df <- cbind(df,skew=skew,kurt=kurt)
##     ## apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
##     ## rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
##     ## data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
## })

## long <- do.call(rbind,by.sim)
## long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
## long <- aggregate(val ~ ., median,data=long)
## by.stats <- split(long,long$stat)
## ## sapply(stats,function(df)coef(lm(log(abs(df$val))~log(df$n)))[2])
## by.stats <- lapply(by.stats, function(df)
##     rates <- sapply(split(df,list(df$skew,df$kurt)), function(smalldf) {
##     if(nrow(smalldf)==0)return(c(NA,NA,NA))
##     ## rate <- if(nrow(smalldf)!=0) {
##     rate <- with(smalldf, coef(lm(log(abs(val))~log(n)))[[2]])
##     ## } else {
##         ## NA
##     ## }
##     c(skew=smalldf$skew[1],kurt=smalldf$kurt[1],rate=rate)
## }) |> t()
## )
## df <- rbind(data.frame(stat='ks.LR',by.stats[['ks.LR']]), data.frame(stat='ks.LR.bartlett',by.stats[['ks.LR.bartlett']]))
## rownames(df) <- NULL
## df <- na.omit(df)

## library(lattice)
## levelplot(rate ~ skew * kurt | stat, data = df,
##           contour = TRUE,      
##           labels = TRUE,       
##           pretty = TRUE,       
##           region = TRUE,       
##           col.regions = gray.colors(100, start = 0.9, end = 0.1),
##           layout = c(2, 1),    
##           xlab = "Skewness", 
##           ylab = "Excess Kurtosis",
##           panel = function(...) {
##               panel.levelplot(...)
##               panel.points(0, 0, pch = 4, col = "red", cex = 1.5) 
##           })

## subdf <- subset(df,stat=='ks.LR.bartlett' & abs(skew)<=.4 & abs(kurt)<=.4)
## matrix(subdf$rate,nrow=length(unique(subdf$skew)))






## ## skew and kurt separately instead of cartesian product

## ## process  rate sim (IID data, using R with gauss-hermite)
## ## save.path <- './sim/robustness/'
## save.path <- './sim/robustness/poisson/'
## save.path <- './sim/robustness/binomial-irr/'
## filelist <- dir(save.path)
## filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
## print(length(filelist))
## by.sim <- lapply(filelist, function(file) {
##     load(paste0(save.path,file))
##     df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
##     df <- do.call(rbind,df.lst)
##     df <- cbind(df,skew=skew,kurt=kurt)
##     ## apply(by.n,1,function(ks)coef(lm(log(ks)~log(ns)))[2])
##     ## rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
##     ## data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
## })
## long <- do.call(rbind,by.sim)
## long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
## long <- aggregate(val ~ ., median,data=long)
## by.stats <- split(long,long$stat)
## ## sapply(stats,function(df)coef(lm(log(abs(df$val))~log(df$n)))[2])
## by.stats <- lapply(by.stats, function(df) 
##     ## rates <- sapply(split(df,list(df$skew,df$kurt)), function(smalldf)
##     rates <- sapply(split(df,paste0(df$skew,df$kurt)), function(smalldf) {
##         if(nrow(smalldf)==0)return(c(NA,NA,NA))
##         ## rate <- if(nrow(smalldf)!=0) {
##         rate <- with(smalldf, coef(lm(log(abs(val))~log(n)))[[2]])
##         ## } else {
##         ## NA
##         ## }
##         c(skew=smalldf$skew[1],kurt=smalldf$kurt[1],rate=rate)
##     }) |> t()    
##     )
## df <- rbind(data.frame(stat='ks.LR',by.stats[['ks.LR']]), data.frame(stat='ks.LR.bartlett',by.stats[['ks.LR.bartlett']]))
## rownames(df) <- NULL
## ## df <- na.omit(df)

## df <- df[order(df$skew,df$kurt),]
## with(subset(df,kurt==0 & stat=='ks.LR'),  plot(rate ~ skew,type='l',ylim=range(df$rate)))
## with(subset(df,kurt==0 & stat=='ks.LR.bartlett'),  lines(rate ~ skew,col=2))

## df <- df[order(df$skew,df$kurt),]
## with(subset(df,skew==0 & stat=='ks.LR'),  plot(rate ~ kurt,type='l',ylim=range(df$rate)))
## with(subset(df,skew==0 & stat=='ks.LR.bartlett'),  lines(rate ~ kurt,col=2))




## compute rates, then aggregate
## process  rate sim (IID data, using R with gauss-hermite)
save.path <- './sim/robustness/binomial-prop2/'
## save.path <- './sim/robustness/poisson/'
## save.path <- './sim/robustness/binomial-irr/'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    load(paste0(save.path,file))
    ## df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
    ## df <- do.call(rbind,df.lst)
    ## df <- cbind(df,skew=skew,kurt=kurt)
    rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
    data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
})
long <- do.call(rbind,by.sim)
long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
rownames(long) <- NULL
long.agg <- aggregate(rate ~ ., median,data=long)
long.agg <- subset(long.agg,subset=stat=='ks.LR' | stat=='ks.LR.bartlett')
plot(rate~skew, data=subset(long.agg,stat=='ks.LR'&kurt==0),type='l',ylim=range(long.agg$rate))
lines(rate~skew, data=subset(long.agg,stat=='ks.LR.bartlett'&kurt==0),col=2)



#########################################################333


## varying tau2 also
## compute rates, then aggregate
## process  rate sim (IID data, using R with gauss-hermite)
## save.path <- './sim/robustness/binomial-prop/'
## save.path <- './sim/robustness/poisson/'
## save.path <- './sim/robustness/binomial-irr/'
## save.path <- './sim/robustness/binomial-irr-tau2/'
## save.path <- './sim/robustness/hyperg/'
## save.path <- './sim/260819/'
save.path <- './sim/bn-prop/'
## save.path <- './sim/bn-irr/'
## save.path <- './sim/hn/'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    ## print(file)   
    load(paste0(save.path,file))
    ## if(is.character(by.n[1]))bad <<- c(bad,file)
    ## print(by.n)
    ## df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
    ## df <- do.call(rbind,df.lst)
    ## df <- cbind(df,skew=skew,kurt=kurt)
    rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
    data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt,tau2.null=tau2.null) 
    ## data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
})
long <- do.call(rbind,by.sim)
long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
rownames(long) <- NULL
long.agg <- aggregate(rate ~ ., mean,data=long)
long.agg <- subset(long.agg,subset=stat=='ks.LR' | stat=='ks.LR.bartlett')
## long.agg <- subset(long.agg,tau2.null%in%c(.6889))
df <- subset(long.agg,stat=='ks.LR.bartlett')
tau2.nulls <- unique(df$tau2.null)
col <- gray.colors(length(tau2.nulls)) |> rev()
## png('./figs/sim-tau2.png')
plot(NA,type='n',xlim=range(df$skew),ylim=c(-2.1,max(long.agg$rate)),xlab='skew',ylab='log rate')
for(i in seq_along(tau2.nulls))
    lines(rate ~ skew, data=subset(df,tau2.null==tau2.nulls[i]),col=col[i])
df <- subset(long.agg,stat=='ks.LR') 
tau2.nulls <- unique(df$tau2.null)
col <- gray.colors(length(tau2.nulls)) |> rev()
## plot(NA,type='n',xlim=range(df$skew),ylim=c(-2,max(df$rate)))
for(i in seq_along(tau2.nulls))
    lines(rate ~ skew, data=subset(df,tau2.null==tau2.nulls[i]),col=col[i],lty=2)
legend('bottomright',col=col,lty=1,legend=sapply(tau2.nulls,\(tau2.null)bquote(tau^2==.(tau2.null))))
## dev.off()



## generate table

rate.sims<- lapply( c(bn.prop='./sim/bn-prop/',
                      ## hn='./sim/robustness/binomial-prop2/',## !!!placeholder
                      hn='./sim/hn/',## !!!placeholder
                      pn='./sim/robustness/poisson/',
                      bn.irr='./sim/bn-irr/'), 
          function(save.path) {
              filelist <- dir(save.path)
              filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
              print(length(filelist))
              by.sim <- lapply(filelist, function(file) {
                  load(paste0(save.path,file))
                  ## df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
                  ## df <- do.call(rbind,df.lst)
                  ## df <- cbind(df,skew=skew,kurt=kurt)
                  rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
                  ## data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt,tau2.null=tau2.null)
                  data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
              })
              long <- do.call(rbind,by.sim)
              long$skew <- round(long$skew,2);long$kurt <- round(long$kurt,2)
              long <- subset(long,kurt==0L)
              long$kurt <- NULL
              rownames(long) <- NULL
              ## long.agg <- aggregate(rate ~ ., median,data=long)
              long.agg <- aggregate(rate ~ ., data=long, FUN=function(x)c(median=median(x),sd=sd(x)))
              long.agg <- subset(long.agg,subset=stat=='ks.LR' | stat=='ks.LR.bartlett')
          })
rate.sims <- lapply(rate.sims,function(df)do.call(data.frame,df))
## temp filler vals!!!!!!!!!
## rate.sims['bn.irr'] <-
    ## rate.sims['hyperg'] <- rate.sims[2]
## rate.sims$bn.prop <- rbind(rate.sims$bn.prop[1:2,],rate.sims$bn.prop)
## rate.sims$bn.prop[1:2,'skew'] <- 0.2
rate.sims <- lapply(rate.sims,function(df)df[df$skew>=0,])

## # 1. Combine median and SD
## df$val <- sprintf("%.2f (%.2f)", df$rate.median, df$rate.sd)
## # 2. Pivot
## wide_df <- reshape(
##   df[, c("stat", "skew", "val")], 
##   idvar = "stat", 
##   timevar = "skew", 
##   direction = "wide"
## )
## ## Clean up headers
## skews <- unique(df$skew)
## wide_df$stat <- ifelse(wide_df$stat == "ks.LR", "LR (standard)", "LR (Bartlett)")
## ## 3. Print directly as a LaTeX table string
## sink('./figs/rate-table.tex')
## cat("\\begin{table}[htbp]\n\\centering\n\\begin{tabular}{l", rep("c", length(skews)), "}\n\\toprule\n", sep = "")
## cat("Statistic & ", paste(skews, collapse = " & "), " \\\\\n\\midrule\n", sep = "")
## for(i in 1:nrow(wide_df)) {
##   vals <- unlist(wide_df[i, -1])
##   cat(wide_df$stat[i], " & ", paste(vals, collapse = " & "), " \\\\\n", sep = "")
## }
## cat("\\bottomrule\n\\end{tabular}\n\\end{table}\n")
## sink()

generate.stacked.latex <- function(study.dfs) {
    ## Helper to process a single study data frame
    process.study <- function(study.name, df) {
        ## Failsafe: if 'rate' is a matrix (from aggregate), flatten it first
        if ("rate" %in% names(df) && is.matrix(df$rate)) {
            df <- do.call(data.frame, df)
        }
        ## Drop skew = 0.25
        df <- df[df$skew != 0.25, ]
        ## Format cell values
        df$val <- sprintf("%.2f (%.2f)", df$rate.median, df$rate.sd)
        ## Pivot wide
        wide <- reshape(
            df[, c("stat", "skew", "val")], 
            idvar = "stat", 
            timevar = "skew", 
            direction = "wide"
        )
        ## Format statistic names
        wide$stat <- ifelse(wide$stat == "ks.LR", "LR (standard)", "LR (Bartlett)")
        ## Add study name column
        study.name <- switch(study.name,
            bn.prop = '\\eqref{mod:bn-prop}',
            hn = '\\eqref{mod:hn}',
            pn = '\\eqref{mod:pn}',
            bn.irr = '\\eqref{mod:bn-irr}',
            )
        wide$study <- study.name
        return(wide)
    }
    ## Process and stack all data frames
    processed.list <- lapply(names(study.dfs), function(nm) process.study(nm, study.dfs[[nm]]))
    stacked.df <- do.call(rbind, processed.list)
    ## Extract skews from the wide column names (val.0, val.0.05, etc.)
    val.cols <- grep("^val\\.", names(stacked.df), value = TRUE)
    ## Strip "val." to get numeric skews and order them properly
    skews.numeric <- as.numeric(gsub("^val\\.", "", val.cols))
    sort.order <- order(skews.numeric)
    val.cols <- val.cols[sort.order]
    skews <- skews.numeric[sort.order]
    ## Total number of columns = Study + Statistic + length(skews)
    total.cols <- 2 + length(skews)
    ## Print LaTeX Table Output
    ## cat("\\begin{table}[htbp]\n\\centering\n")
    cat("\\begin{tabular}{l l", rep("c", length(skews)), "}\n\\toprule\n", sep = "")
    cat("& & \\multicolumn{5}{c}{\\textbf{Skew of random effect}} \\\\")
    cat("\\cmidrule(lr){3-7}")
    cat("\\textbf{Model} & \\textbf{Statistic} & ", paste(skews, collapse = " & "), " \\\\\n\\midrule\n", sep = "")
    ## Loop over studies to format rows
    study.names <- unique(stacked.df$study)
    for (i in seq_along(study.names)) {
        s.name <- study.names[i]
        sub.df <- stacked.df[stacked.df$study == s.name, ]
        n.rows <- nrow(sub.df)
        for (r in 1:n.rows) {
            vals <- unlist(sub.df[r, val.cols])
            ## Replace any NAs with a dash just in case a study lacked a specific skew
            vals[is.na(vals)] <- "-"
            row.str <- paste(vals, collapse = " & ")
            if (r == 1) {
                cat(sprintf("\\multirow{%d}{*}{%s} & %s & %s \\\\\n", n.rows, s.name, sub.df$stat[r], row.str))
            } else {
                cat(sprintf(" & %s & %s \\\\\n", sub.df$stat[r], row.str))
            }
        }
        ## Dynamic horizontal divider width
        if (i < length(study.names)) {
            cat(sprintf("\\cmidrule(lr){1-%d}\n", total.cols))
        }
    }
    cat("\\bottomrule\n\\end{tabular}\n")
    ## cat("\\caption{Simulation results. The cells contain the medians and SDs of the log of the rate of the empirical Kolmogorov-Smirnov distance to a chi-squared distribution.}\n")
    ## cat("\\label{tab:simulation}\n")
    ## cat("\\end{table}\n")
}
sink('./figs/rate-table.tex')
generate.stacked.latex(rate.sims)
sink()




## varying tau2, with loess


save.path <- './sim/bn-irr/'
## save.path <- './sim/hn/'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    ## print(file)   
    load(paste0(save.path,file))
    rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
    data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt,tau2.null=tau2.null) 
    ## data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
})

long <- do.call(rbind,by.sim)
long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
rownames(long) <- NULL
long.agg <- aggregate(rate ~ ., median,data=long)
long.agg <- subset(long.agg,subset=stat=='ks.LR' | stat=='ks.LR.bartlett')
## long.agg <- subset(long.agg,tau2.null%in%c(.1,.2,.3))
df <- subset(long.agg,stat=='ks.LR.bartlett')
tau2.nulls <- unique(df$tau2.null)
col <- gray.colors(length(tau2.nulls)) |> rev()
## png('./figs/sim-tau2.png')
plot(NA,type='n',xlim=range(df$skew),ylim=c(-2.1,max(long.agg$rate)),xlab='skew of random effect',ylab='log rate')
for(i in seq_along(tau2.nulls))
    ## lines(rate ~ skew, data=subset(df,tau2.null==tau2.nulls[i]),col=col[i])
    with(subset(df,tau2.null==tau2.nulls[i]), lines(skew, predict(loess(rate ~ skew)), col=col[i]))
df <- subset(long.agg,stat=='ks.LR')
tau2.nulls <- unique(df$tau2.null)
col <- gray.colors(length(tau2.nulls)) |> rev()
## plot(NA,type='n',xlim=range(df$skew),ylim=c(-2,max(df$rate)))
for(i in seq_along(tau2.nulls))
    ## lines(rate ~ skew, data=subset(df,tau2.null==tau2.nulls[i]),col=col[i],lty=2)
    with(subset(df,tau2.null==tau2.nulls[i]), lines(skew, predict(loess(rate ~ skew)), col=col[i],lty=2))
legend('bottomright',col=col,lty=1,legend=sapply(tau2.nulls,\(tau2.null)bquote(tau^2==.(tau2.null))))
## dev.off()









## further anlaysis for hn and bn-irr models

save.path <- './sim/bn-irr/'
## save.path <- './sim/hn/'
filelist <- dir(save.path)
filelist <- filelist[grep('^save[-0-9]+\\.RData',filelist)]
print(length(filelist))
by.sim <- lapply(filelist, function(file) {
    ## print(file)   
    load(paste0(save.path,file))
    ## if(is.character(by.n[1]))bad <<- c(bad,file)
    ## print(by.n)
    ## df.lst <- lapply(1:nrow(by.n),function(i)data.frame(stat=rownames(by.n)[i],val=by.n[i,],n=ns))
    ## df <- do.call(rbind,df.lst)
    ## df <- cbind(df,skew=skew,kurt=kurt)
    rates <- apply(by.n,1,function(x)coef(lm(log(abs(x))~log(ns)))[2])
    data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt,tau2.null=tau2.null) 
    ## data.frame(rate=rates,stat=names(rates),skew=skew,kurt=kurt)
})
long <- do.call(rbind,by.sim)
long$skew <- round(long$skew,2); long$kurt <- round(long$kurt,2)
rownames(long) <- NULL
df <- long.agg <- aggregate(rate ~ ., data=long, FUN=function(x)c(median=median(x),sd=sd(x)))
df <- cbind(df[,names(df)!='rate'],df$rate)
df <- df[df$skew>=0,]

generate_latex_tabular <- function(df, digits = 2, tau_digits = 1) {
  # 1. Filter dataframe to include only the required statistics
  target_stats <- c("ks.LR", "ks.LR.bartlett")
  df <- df[df$stat %in% target_stats, ]
  
  # 2. Map statistic names to LaTeX labels
  stat_map <- c(
    "ks.LR"          = "LR (original)",
    "ks.LR.bartlett" = "LR (Bartlett)"
  )
  
  # 3. Format median and SD into "median (sd)" strings
  fmt <- paste0("%.", digits, "f (%.", digits, "f)")
  df$cell_val <- sprintf(fmt, df$median, df$sd)
  
  # 4. Extract unique factor levels maintaining grid ordering
  tau_vals  <- unique(df$tau2.null)
  stat_vals <- intersect(target_stats, unique(df$stat))
  skew_vals <- sort(unique(df$skew))
  
  num_skew  <- length(skew_vals)
  skew_cols <- paste(skew_vals, collapse = " & ")
  
  # 5. Construct LaTeX tabular environment
  lines <- c()
  lines <- c(lines, sprintf("\\begin{tabular}{ll%s}", paste(rep("c", num_skew), collapse = "")))
  lines <- c(lines, "\\toprule")
  lines <- c(lines, sprintf(" &  & \\multicolumn{%d}{c}{\\textbf{Skew of random effect}} \\\\", num_skew))
  lines <- c(lines, sprintf("\\cmidrule(lr){3-%d}", num_skew + 2))
  lines <- c(lines, sprintf("\\textbf{$\\tau^2$} & \\textbf{Statistic} & %s \\\\", skew_cols))
  lines <- c(lines, "\\midrule")
  
  # Loop through tau2 levels and statistic types
  for (i in seq_along(tau_vals)) {
    t_val <- tau_vals[i]
    
    # Format tau2 to 1 decimal place
    formatted_tau <- sprintf(paste0("%.", tau_digits, "f"), t_val)
    ## tau_label <- sprintf("$\\tau^2 = %s$", formatted_tau)
    tau_label <- sprintf("$%s$", formatted_tau)
    
    for (j in seq_along(stat_vals)) {
      s_val <- stat_vals[j]
      s_label <- if (s_val %in% names(stat_map)) stat_map[[s_val]] else s_val
      
      # Gather cell values across skew levels
      row_cells <- sapply(skew_vals, function(sk) {
        sub <- df[df$tau2.null == t_val & df$stat == s_val & df$skew == sk, ]
        if (nrow(sub) > 0) sub$cell_val[1] else ""
      })
      
      # Use \multirow on the first row of each tau2 block
      first_col <- if (j == 1) sprintf("\\multirow{%d}{*}{%s}", length(stat_vals), tau_label) else ""
      row_str <- paste(c(first_col, s_label, row_cells), collapse = " & ")
      lines <- c(lines, paste0(row_str, " \\\\"))
    }
    
    if (i < length(tau_vals)) {
      lines <- c(lines, "\\midrule")
    }
  }
  
  lines <- c(lines, "\\bottomrule")
  lines <- c(lines, "\\end{tabular}")
  
  return(paste(lines, collapse = "\n"))
}

latex_code <- generate_latex_tabular(df,tau_digits=2)
cat(latex_code)
sink('./figs/rate-table-bn-irr.tex')
cat(latex_code)
sink()
