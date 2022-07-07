
library(data.table)

# generate a table of events independent across samples
event.table <- function(N,pe) {
# N is the number of samples to simulate
# pe is a vector of probabilities of an event per sample 
  N.events <- length(pe)
  if (is.null(names(pe))) {
    names(pe) <- paste('event_',1:length(pe),sep='')
  }
  et <- data.table(id=paste('sample',1:N,sep=''))
  for (nm in names(pe)) {
    #et[[paste('event',n,sep='')]] <- runif(N) < pe[i]
    et[[nm]] <- runif(N) < pe[[nm]]
  }
  return(et)
}

# perform fisher exact tests on all pairs of events
pairwise.fish <- function(et) {
  et <- et[,-1] # remove id column
  nv = length(et)
  tl = nv*(nv-1)/2
  p.table <- data.table(i=rep(0,tl),j=rep(0,tl))
  k <- 1
  for (i in 1:(nv-1)) {
    for (j in (i+1):nv) {
      p.table$i[k] <- i
      p.table$j[k] <- j
      p.table$p.corr[k] <- fisher.test(et[[i]],et[[j]],alternative='less')$p.value
      p.table$p.anti[k] <- fisher.test(et[[i]],et[[j]],alternative='greater')$p.value
      p.table$p.both[k] <- fisher.test(et[[i]],et[[j]],alternative='two.sided')$p.value
      k <- k + 1
    }
  }
  
  return (p.table)
}

# perform fisher exact test on all pairs of events
pairwise.p.table <- function(et,tail='less') {
  et <- et[,-1]
  #et <- event.table(N,pe)[,-1]
  nv = length(et)
  tl = nv*(nv-1)/2
  p.table <- data.table(i=rep(0,tl),j=rep(0,tl),p=rep(0,tl))
  k <- 1
  for (i in 1:(nv-1)) {
    for (j in (i+1):nv) {
      #!print(paste(i,j,k))
      p.table$i[k] <- i
      p.table$j[k] <- j
      p.table$p[k] <- fisher.test(et[[i]],et[[j]],alternative=tail)$p.value
      k <- k + 1
    }
  }
  
  return (p.table[order(p)])
}

event.p.table <- function(N,pe,tail='less') {
  et <- event.table(N,pe)[,-1]
  return(pairwise.p.table(et,tail=tail))
}

#### Main Code


# probability of event per sample vector
pe <- seq(0.30,0.70,0.02)


sim.nothing <- FALSE

# either simulate different rates of same propensities
# or different propensity patterns
sim.rate.var <- TRUE

if (sim.nothing) {
  et <- event.table(6000,pe)
} else {
  if (sim.rate.var) {
    # simulate samples w/different (proportional) rates of events
    #et <- event.table(6000,pe)
    et1 <- event.table(2000,pe)
    et2 <- event.table(2000,1.1*pe)
    et3 <- event.table(2000,0.9*pe)
    et <- rbind(et1,et2,et3)
  } else {
    # simulate mixed classes of samples with different individual event rates
    variant.fraction = 0.2
    pe2 <- pe
    ix <- sample(seq_along(pe),size=round(variant.fraction*length(pe)))
    pe2[ix] <- sample(pe[ix])
    et <- rbind(event.table(3000,pe),event.table(3000,pe2))
  }
}
# et is the simulated event table we will analyze

# the output tables pt.corr and pt.anti have one-sided fisher P values
# for correlation and anti-correlation, respectively
library(stats)



pt.corr <- pairwise.p.table(et,tail='greater')
pt.corr[,fdr:=p.adjust(p,method='BH')]
pt.anti <- pairwise.p.table(et,tail='less')
pt.anti[,fdr:=p.adjust(p,method='BH')]

pp <- pairwise.fish(et)

### plot the p-values
library(ggplot2)

graygreen <- '#669966'
grayblue <- '#6666CC'
reddish <- '#CC6666'

xxx <- data.table(fisher.p=pt.corr$p,unif.p=(1:nrow(pt.corr))/nrow(pt.corr))

# plot uniform vs observed "KS plot"

alpha <- 0.25 # acceptable false discovery rate
gp <- ggplot(xxx,aes(x=fisher.p,y=unif.p)) + geom_line(color=grayblue) + 
    geom_abline(slope=1,intercept=0,color='gray') + 
    geom_rug(sides='b',color=grayblue) +
    geom_abline(intercept=0,slope=1/alpha,color=graygreen) +
    geom_abline(intercept=1-1/alpha,slope=1/alpha,color=graygreen) +
    geom_text(aes(label=paste('FDR < ',alpha,sep=''),x=0.5*alpha,y=0.5),
              nudge_x=-0.02,angle=atan(1/alpha)*180/pi,color=graygreen)

#pdf(file='~/R/ks.pdf')
print(gp)
#dev.off()

# draw a classic q-q plot (correlation only)

pt.both <- rbind(pt.corr[,side := 'corr'],pt.anti[,side:='anti'] )
pt.both[,log10.fisher.p := -log(p,base=10)]
pt.both <- pt.both[order(p),]
pt.both[,log10.unif.p := -log((1:.N)/.N,base=10)]
pt.both[,fdr := p.adjust(p,method='BH')]

# add conf intervals
pt.both[,k := 1:.N]
pt.both[,p0 := -log(qbeta(0.05,k+1,.N-k+1),base=10)]
pt.both[,p1 := -log(qbeta(0.95,k+1,.N-k+1),base=10)]

#!corr.p <- data.table(log10.fisher.p=-log(pt.corr$p,base=10),log10.unif.p=-log(1:nrow(pt.corr)/nrow(pt.corr),base=10),side='corr')
#!mx <- ceiling(max(-log(c(pt.corr$p,1/nrow(pt.corr)),base=10)))
#!anti.p <- data.table(log10.fisher.p=-log(pt.anti$p,base=10),log10.unif.p=-log(1:nrow(pt.anti)/nrow(pt.anti),base=10),side='anti')
#!mx <- ceiling(max(rbind(corr.p,anti.p)$log10.fisher.p))

#gp1 <- ggplot(rbind(corr.p,anti.p),aes(x=log10.unif.p,y=log10.fisher.p,color=side)) + geom_point(alpha=0.5) +
  gp1 <- ggplot(pt.both,aes(x=log10.unif.p,y=log10.fisher.p,color=side)) + geom_point(alpha=0.7) +
  geom_abline(slope=1,intercept=0,color='gray') + expand_limits(x=c(0,mx),y=c(0,mx)) +
    geom_abline(slope=1,intercept=-log(alpha,base=10),color=graygreen) +
    geom_text(aes(label=paste('FDR < ',alpha,sep=''),x=-log(alpha,base=10),y=2*-log(alpha,base=10)),
              nudge_y=0.04,nudge_x=-0.04,angle=45,color=graygreen) +
    geom_ribbon(aes(ymin=p1,ymax=p0),color=graygreen,alpha=0.2)
    
#pdf(file='~/R/qq.pdf')
print(gp1)
#dev.off()
