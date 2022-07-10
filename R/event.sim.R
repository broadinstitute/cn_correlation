
library(data.table)

# generate a table of events independent across samples
event.table <- function(N,pe) {
  # process vector arguments an element at a time
  if (length(N) > 1) {
    if (length(N) != length(pe)) stop('vectors must be of equal length')
    et <- rbindlist(lapply(seq_along(N), function(i) event.table(N[i],pe[[i]])))
  }
  else {
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
  }
  return(et)
}

# perform fisher exact tests on all pairs of events in an event table
# in: et - event table
#     side - which test(s) 'corr', 'anti', 'both', or a list of any combination of these
# 
pairwise.fisher.test <- function(et,side='both') {
  et <- et[,-1] # remove simple id column
  # build pairwise table
  nv = length(et) # number of events
  pairx <- which(upper.tri(matrix(1:(nv*nv),nrow=nv)),arr.ind=TRUE)
  p.table <- data.table(i=pairx[,'row'],j=pairx[,'col'])
  for (k in 1:nrow(p.table)) {
    if (is.element('corr',side)) {
      p.table$corr[k] <- fisher.test(et[[p.table$i[k]]],et[[p.table$j[k]]],alternative = 'greater')$p.value
    }
    if (is.element('anti',side)) {
      p.table$anti[k] <- fisher.test(et[[p.table$i[k]]],et[[p.table$j[k]]],alternative = 'less')$p.value
    }
    if (is.element('both',side)) {
      p.table$both[k] <- fisher.test(et[[p.table$i[k]]],et[[p.table$j[k]]],alternative = 'two.sided')$p.value
    }
  }
  return (p.table)
}


#' Swaps required to adjust odds ration of a contingency table 
#' @description 
#' Based on marginals of the supplied contingency table of event associations, 
#' calculate the number of marginal-conserving swaps required to adjust the odds
#' ratio of a table of independent event associations with the same marginals.
#' 
#' @param t1 2x2 contingency table of co-occurrences
#' @param rho desired odds ratio, > 1 for correlation, < 1 for anti-correlation
#' @return Number of swaps, > 0 for correlation, < 0 for anti-correlation.
#'  The return value is in general a non-integer, and will be infinite if the 
#'  odds ration cannot be achieved by swapping.
#' 
swaps4odds <- function(t1,rho) {
  # create independent table (a,b,c,d) from marginals
  tm <- addmargins(t1)
  N <- tm['Sum','Sum']
  a <- tm['FALSE','Sum'] * tm['Sum','FALSE'] / N
  b <- tm['FALSE','Sum'] * tm['Sum','TRUE']  / N
  c <- tm['TRUE','Sum']  * tm['Sum','FALSE'] / N
  d <- tm['TRUE','Sum']  * tm['Sum','TRUE']  / N
  
  # find x that will achieve the desired odds ratio
  x <- polyroot( c(a*d-rho*b*c, a+d+rho*c+rho*d, (1-rho)))
  x <- x[abs(x) == min(abs(x))]
  if (length(x) > 1) {
    x <- x[1]
  }
  
  # complex roots indicate odds ratio cannot be achieved
  if (!isTRUE(all.equal(Im(x),0))) {
    # return correctly signed infinity 
    x <- Inf * sign(log(rho))
  }
  # return root with smaller absolute value
  return(Re(x))
}

#' Create correlations/anti-correlations in event table
#' @param et data.table of events
#' @param corr.tab data.table specifying correlations to create
#' for selected pairs of events
#' @return data.table of adjusted events with same marginals
create.corr <- function(et,corr.tab) {
  id <- et[,1]
  et <- et[,-1]
  et0 <- et
  for (i in seq_along(corr.tab$e1)) {
    swap.cols <- c(corr.tab$e1[i],corr.tab$e2[i])
    pair.vals <- as.matrix(et[,..swap.cols])
    ### simulate correlation
    if (corr.tab$log2.fold[i] > 0) {
      # find samples with both types of anti-correlated event pairs
      bs <- which(pair.vals[,1] & !pair.vals[,2])
      cs <- which(!pair.vals[,1] & pair.vals[,2])
      # randomly shorten the longer type to same size
      if (length(bs)>length(cs)) {
        bs <- sample(bs,length(cs))
      }
      else {
        cs <- sample(cs,length(bs))
      }
      # calculate the number of swaps needed to achieve target odds ratio given marginals
      swaps.needed <- swaps4odds(table(et[,..swap.cols]),2^corr.tab$log2.fold[i])
      # convert to a probability such that expected number of swaps is target
      if (is.infinite(swaps.needed)) {
        p.swap <- 1.0;
      }
      else {
        p.swap <- min(1.0, swaps.needed / length(bs))
      #!print(paste(corr.tab$e1[i],corr.tab$e2[i],p.swap))
      }
      if (p.swap == 1.0) {
        warning(paste('swaps saturated columns ',swap.cols[1],' x ',swap.cols[2],'\n',sep=''))
      }
      for (j in seq_along(bs)) {
        if (runif(1) < p.swap) {
          # for correlation, swap event choices are where event statuses differbetween samples, excepting the swap.cols
          swap.choices <- which( (et[bs[j],]!=et[cs[j],]) & !is.element(colnames(et),swap.cols) )
          if (length(swap.choices) > 0) {
            swap.x <- sample(swap.choices,1)
            if (et[bs[j],..swap.x] == TRUE) {
              # swap outside column with column 2 of pair
              et[bs[j],swap.x] <- FALSE
              et[bs[j],swap.cols[2]] <- TRUE
              et[cs[j],swap.cols[2]] <- FALSE
              et[cs[j],swap.x] <- TRUE
            }
            else {
              # swap outside column with column 1 of pair
              et[bs[j],swap.x] <- TRUE
              et[bs[j],swap.cols[1]] <- FALSE
              et[cs[j],swap.cols[1]] <- TRUE
              et[cs[j],swap.x] <- FALSE
            }
          }
          else {
            warning('no swap choices')
          }
        }
      }
    }
    ### simulate anti-correlation
    if (corr.tab$log2.fold[i] < 0) {
      # find both types of correlated event pairs
      as <- which(!pair.vals[,1] & !pair.vals[,2]) ##!
      ds <- which(pair.vals[,1] & pair.vals[,2]) ##!
      # randomly shorten the longer type to same size
      if (length(ds)>length(as)) { ##
        ds <- sample(ds,length(as)) ##
      }
      else {
        as <- sample(as,length(ds)) ##
      }
      # calculate the number of swaps needed to achieve target odds ratio given marginals
      swaps.needed <- - swaps4odds(table(et[,..swap.cols]),2^corr.tab$log2.fold[i]) ##! sign
      # convert to a probability such that expected number of swaps is target
      if (is.infinite(swaps.needed)) {
        p.swap <- 1.0;
      }
      else {
        p.swap <- min(1.0, swaps.needed / length(as)) ##
        #!print(paste(corr.tab$e1[i],corr.tab$e2[i],p.swap))
      }
      if (p.swap == 1.0) {
        warning(paste('swaps saturated columns',swap.cols))
      }
      for (j in seq_along(as)) {
        if (runif(1) < p.swap) {
          # find events with a 1 for the as and a 0 for the ds
          swap.choices <- which((et[as[j],]==1) & (et[ds[j],]==0))
          # pick one of the columns of interest
          swap.col <- sample(swap.cols,1)
          # perform the swap
          if (length(swap.choices) > 0) {
            et[as[j],swap.x] <- FALSE
            et[as[j],swap.col] <- TRUE
            et[ds[j],swap.col] <- FALSE
            et[cs[j],swap.x] <- TRUE
          }
          else {
            warning('no swap choices')
          }
        }
      }
    }
    
    # put altered columns back in event table
#!    et[,swap.cols[[1]]] <- pair.vals[,1]
#!    et[,swap.cols[[2]]] <- pair.vals[,2]
    # (so much for the data.table's much-touted parallelism)
    #!print(table(et[,..swap.cols])-table(et0[,..swap.cols]))
  } # end corr.tab loop
  return(cbind(id,et))
}

#### Main Code


# independent uniform pseudo-randomly generated events


library(stats)


### plot the p-values
library(ggplot2)

#!corr.p <- data.table(log10.fisher.p=-log(pt.corr$p,base=10),log10.unif.p=-log(1:nrow(pt.corr)/nrow(pt.corr),base=10),side='corr')
#!mx <- ceiling(max(-log(c(pt.corr$p,1/nrow(pt.corr)),base=10)))
#!anti.p <- data.table(log10.fisher.p=-log(pt.anti$p,base=10),log10.unif.p=-log(1:nrow(pt.anti)/nrow(pt.anti),base=10),side='anti')
#!mx <- ceiling(max(rbind(corr.p,anti.p)$log10.fisher.p))

qq.loglog <- function(pp,alpha=0.25,conf=0.95) {
  # strip to columns of interest and melt
  keepcols <- intersect(names(pp),c('i','j','corr','anti','both'))
  p3 <- melt(pp[,..keepcols],id=c('i','j'),variable.name = 'side',value.name = 'p')
  
  p3[side=='corr', k := order(order(p))]
  p3[side=='anti', k := order(order(p))]
  p3[side=='both', k := order(order(p))]
#!  p3[side=='corr', fdr := p.adjust(p,method='BH')]
#!  p3[side=='anti', fdr := p.adjust(p,method='BH')]
#!  p3[side=='both', fdr := p.adjust(p,method='BH')]
  # create log probabilities
  p3[,log10.unif.p := -log(k/nrow(pp),base=10)]
  p3[,log10.fisher.p := -log(p,base=10)]
  
  graygreen <- '#99BB99'
  grayblue <- '#6666CC'
  reddish <- '#CC6666'

  # add confidence intervals
  p3[,cmin := -log(qbeta((1-conf)/2,k+1,nrow(pp)-k+1),base=10)]
  p3[,cmax := -log(qbeta((1+conf)/2,k+1,nrow(pp)-k+1),base=10)]
  mx <- max(p3$log10.fisher.p)
  gp1 <- ggplot(p3,aes(x=log10.unif.p,y=log10.fisher.p,color=side)) + geom_point(alpha=0.7) +
  geom_abline(slope=1,intercept=0,color='gray') + #!expand_limits(x=c(0,mx),y=c(0,mx)) +
    geom_abline(slope=1,intercept=-log(alpha,base=10),color=graygreen) +
    geom_text(aes(label=paste('FDR < ',alpha,sep=''),x=-log(alpha,base=10),y=2*-log(alpha,base=10)),
              nudge_y=0.04,nudge_x=-0.04,angle=45,color=graygreen) +
    geom_ribbon(aes(ymin=cmax,ymax=cmin),color='gray',alpha=0.2)
  return(gp1)
}

# probability of event vector
pe <- seq(0.78,0.98,0.005)


# uniform, 6000 samples

et <- event.table(6000,pe)

# three populations of samples with different rates
# (is scaling the probabilities linearly the right thing to do here?)
# (seems like we should "age" the probabilities towards 1 to be on model)

#!age.p <- function(p,lambda) return(1 - (1-p)^lambda)
#!et <- event.table(c(2000,2000,2000),list(age.p(pe,0.9),pe,age.p(pe,1.1)))

# two different types of event source
# alpha dials up the degree of difference

#!beta <- 0.6
#!pe2 <- beta*sample(pe) + (1-beta)*pe;
#!gp1 <- ggplot(data.table(pe1=pe,pe2=pe2),aes(x=pe1,y=pe2)) + geom_point()
#!print(gp1)
#!et <- event.table(c(3000,3000),list(pe,pe2))

# et is the simulated event table we will analyze



# simulate event table
#!et <- event.table(6000,pe)


# calculate p-values according to fisher exact test
pp <- pairwise.fisher.test(et,c('corr','anti'))

# plot the results
#pdf(file='~/R/qq.pdf')
print(qq.loglog(pp))
#dev.off()

# engineer specific correlations
corr.tab <- data.table(e1=c('event_1','event_11','event_3','event_13'),
                       e2=c('event_2','event_12','event_4','event_14'),
                       log2.fold=c(0.3,0.3,-0.3,-0.3))
et2 <- create.corr(et,corr.tab)

pp2 <- pairwise.fisher.test(et2,c('corr','anti'))
print(qq.loglog(pp2))

pp2[,q.corr := p.adjust(corr,method='BH')]
pp2[,q.anti := p.adjust(anti,method='BH')]

