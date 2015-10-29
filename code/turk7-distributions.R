library(nullabor)
library(ggplot2)
library(plyr)
library(reshape)
library(fpc)
library(tourr)
library(tidyr)
library(dplyr)

##################
# file of all separation indices for all of the turk7 lineups (is re-calculated in the Rnw). 
metrics <- read.table("arxiv submission/largep-metrics.txt", header=TRUE)




#### New Set of 19 nulls

generate_nulls_2d<-function(n=30, p, noise=1, m=19){
  x<-matrix(rnorm(p*n),ncol=p)
#   if(noise==0){
#     x[1:10,(p-1)]<-x[1:10,(p-1)]+3
#     x[11:20,(p-1)]<-x[11:20,(p-1)]-3
#     x[21:30,p]<-x[21:30,p]+sqrt(27)
#   }
  colnames(x)<-paste("X",1:(p),sep="")
  x<-scale(x)
  x<-data.frame(x, cl=factor(c(rep(1,n/3),rep(2,n/3),rep(3,n/3))))
  d=2
  
#   optima <- save_history(x[,-(p+1)], tour_path=guided_tour(index_f=pda_pp(cl=x[,(p+1)], lambda=0.2), max.tries=1000), max_bases=100, rescale=F)
#   nbases<-dim(optima)[3]
#   optima.global<-unclass(optima)[,,nbases]
#   
#   projdata.true<-data.frame(as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)], nbases=rep(nbases,n))
  
  projdata.samples<-NULL
  flag <- 0
  while(flag < m) {
    x[,(p+1)]<-sample(x[,(p+1)])
    optima <- save_history(x[,-(p+1)], guided_tour(index_f=pda_pp(cl=x[,(p+1)], lambda=0.2), max.tries=250), max_bases=200, rescale=F)
    nbases<-dim(optima)[3]
    optima.global<-unclass(optima)[,,nbases]
    projdata<-data.frame(as.matrix(x[,-(p+1)])%*%optima.global, cl=x[,(p+1)], nbases=rep(nbases,30))
    lamb<-summary(manova(cbind(X1, X2)~cl, data=projdata), test="Wilks")[[4]][3]
    if(lamb < 0.001){
      flag <- flag + 1
      projdata$.samples <- flag
      projdata.samples<-rbind(projdata.samples, projdata)
    } else {
      projdata <- NULL
  #    projdata.samples <- rbind(projdata.samples, projdata)
      flag = flag
    }
    cat(flag, lamb, "\n") 
  }
#  projdata.samples$.n <- rep(1:(m-1), each = n)
  #pos<-sample(m,1)
#  lineup.data<-lineup(true=projdata.true, samples=projdata.samples,pos=1)
  return(projdata.samples)
}

n <- 30; p <- 80
for (i in 1:9) {
  nulli <- generate_nulls_2d(p = 80, n = 30, m =19)
  write.csv(nulls2, file=sprintf("nulls-turk7/null-0%d.csv", i), row.names=FALSE)
}

for (i in 10:99) {
  nulli <- generate_nulls_2d(p = 80, n = 30, m =19)
  write.csv(nulls2, file=sprintf("nulls-turk7/null-%d.csv", i), row.names=FALSE)
}

null_dists <- function(lineup) {
  m <- max(lineup$.samples)
  dat.dbn <- dat.dms <- dat.das <- NULL
  pos1 <- pos2 <- NULL
  for (i in 1:19) {
    dbn <- dms <- das <- NULL
    for (j in 1:19) {
      if (i != j) {
        Xi <- subset(lineup, .samples == i)
        Xj <- subset(lineup, .samples == j)
        dbn = c(dbn, bin_dist(Xi, Xj, lineup.dat = lineup, X.bin = 5, Y.bin = 5))
        dms = c(dms, sep_dist(Xi, Xj, clustering = TRUE, nclust = 3))
        das = c(das, sep_dist(Xi, Xj, clustering = TRUE, 
                       nclust = 3, type="average.toother"))
      }
    }
    pos1 <- c(pos1, i)
    pos2 <- c(pos2, j)
    dat.dbn <- c(dat.dbn, mean(dbn))
    dat.dms <- c(dat.dms, mean(dms))
    dat.das <- c(dat.das, mean(das))
  }
  data.frame(pos1, pos2, dbn = dat.dbn, dms = dat.dms, das = dat.das)
}

setwd("nulls-turk7")
dists <-  list.files(pattern=".csv$") %>% 
  lapply(read.csv, stringsAsFactors = FALSE) %>%
  lapply(null_dists) %>% ldply(function(x) x)
setwd("..")

####
# densities

qplot(data = dists, dbn, geom = "density", fill = I("grey80"), colour = I("grey80"), 
      xlab = "Binned (5,5) density", ylab = "")                                                                                                                     max(density(dat.b)$y) + 0.1))

ggsave("bin-dist-largep-6-4-new-1.pdf", height = 5, width = 5.5)        

### Sep dist

dat$cl <- as.numeric(dat$cl)
ddd <- distmet(dat, var = c("X1", "X2", "cl"), 'sep_dist', null_permute("cl"), pos = 1, dist.arg = list(clustering = TRUE, nclust = 3))  


pos <- 1; m = 20

qplot(dat.s, geom = "density", fill = I("grey80"), colour = I("grey80"), 
      xlab = "Permutation distribution", ylab = "") + geom_segment(aes(x = ddd$lineup$mean.dist[ddd$lineup$plotno != 
                                                                                                  pos], xend = ddd$lineup$mean.dist[ddd$lineup$plotno != pos], y = rep(0.01 * min(density(dat.s)$y), 
                                                                                                                                                                       (m - 1)), yend = rep(0.05 * max(density(dat.s)$y), (m - 1))), size = 1, alpha = I(0.7)) + 
  geom_segment(aes(x = ddd$lineup$mean.dist[ddd$lineup$plotno == pos], xend = ddd$lineup$mean.dist[ddd$lineup$plotno == 
                                                                                                     pos], y = 0.01 * min(density(dat.s)$y), yend = 0.1 * max(density(dat.s)$y)), 
               colour = "darkorange", size = 1) + geom_text(data = ddd$lineup, y = -0.03 * max(density(dat.s)$y), 
                                                            size = 2.5, aes(x = mean.dist, label = plotno)) + ylim(c(-0.04 * max(density(dat.s)$y), 
                                                                                                                     max(density(dat.s)$y) + 0.1))


ggsave("sep-dist-largep-new-1.pdf", height = 5, width = 5.5)  

