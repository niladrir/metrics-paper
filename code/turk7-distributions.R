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
metrics <- read.table("data/largep-metrics.txt", header=TRUE)



#### New Set of 19 nulls

generate_nulls_2d<-function(n=30, p, noise=1, m=19){
  x<-matrix(rnorm(p*n),ncol=p)
   if(noise==0){
     x[1:10,(p-1)]<-x[1:10,(p-1)]+3
     x[11:20,(p-1)]<-x[11:20,(p-1)]-3
     x[21:30,p]<-x[21:30,p]+sqrt(27)
   }
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

n <- 30; p <- 100
for (i in 0:9) {
  nulli <- generate_nulls_2d(p = 100, n = 30, m =19)
  write.csv(nulli, file=sprintf("nulls-turk7-2/null-0%d.csv", i), row.names=FALSE)
}

for (i in 10:99) {
  nulli <- generate_nulls_2d(p = 100, n = 30, m =19)
  write.csv(nulli, file=sprintf("nulls-turk7-2/null-%d.csv", i), row.names=FALSE)
}

null_dists <- function(lineup) {
  m <- max(lineup$.samples)
  dat.dbn <- dat.dms <- dat.das <- dat.dmin <- dat.ddunn <- NULL
  pos1 <- pos2 <- NULL
  for (i in 1:19) {
    dbn <- dms <- das <- ddunn <- dmin <- NULL
    for (j in 1:19) {
      if (i != j) {
        Xi <- subset(lineup, .samples == i)
        Xj <- subset(lineup, .samples == j)
        dbn = c(dbn, bin_dist(Xi, Xj, lineup.dat = lineup, X.bin = 5, Y.bin = 5))
        dmin = c(dmin, sep_dist(Xi, Xj, clustering = TRUE, type = "min.separation"))
        dms = c(dms, sep_dist(Xi, Xj, clustering = TRUE))
        ddunn = c(ddunn, sep_dist(Xi, Xj, clustering = TRUE, type = "dunn"))
        das = c(das, sep_dist(Xi, Xj, clustering = TRUE, 
                       nclust = 3, type="average.toother"))
      }
    }
    pos1 <- c(pos1, i)
#    pos2 <- c(pos2, j)
    dat.dbn <- c(dat.dbn, mean(dbn))
    dat.dms <- c(dat.dms, mean(dms))
    dat.das <- c(dat.das, mean(das))
    dat.dmin <- c(dat.dmin, mean(dmin))
    dat.ddunn <- c(dat.ddunn, mean(ddunn))
  }
#  data.frame(dbn = mean(dat.dbn), dms = mean(dat.dms), das = mean(dat.das))
  pick <- sample(1:19, 1)
#  data.frame(dbn = dat.dbn[pick], dms = dat.dms[pick], das = dat.das[pick])
  data.frame(dbn = dat.dbn, dms = dat.dms, das = dat.das, ddunn = dat.ddunn, dmin = dat.dmin)
}

setwd("nulls-turk7-100-30")
dists <-  list.files(pattern=".csv$") %>% 
  lapply(read.csv, stringsAsFactors = FALSE) %>%
  lapply(null_dists) %>% ldply(function(x) x)
setwd("..")
write.csv(dists, "data/reference-distances-turk7-100-30.csv", row.names=FALSE)

pos <- 20
submetrics <- subset(metrics, pic_name == "plot_large_p_small_n_30_100_0_2_3.png")
submetrics <- subset(submetrics, pos.1 != pos.2 & pos.2 != pos)
dd3 <- submetrics %>% group_by(pos.1) %>% summarize(
  dbn.mean = mean(dbn),
  das.mean =  mean(das),
  dms.mean = mean(das),
  dmin.mean = mean(dmin),
  ddunn.mean = mean(ddunn),
  len=n())

#### plots
# lineup

theme_lineup <- theme(axis.text = element_blank(), 
                      axis.title = element_blank(),
                      axis.ticks = element_blank(),
                      plot.margin=unit(c(0,0,0,0), unit="cm"))


lpt7 <- read.table("data/turk7/dat_large_p_small_n_30_100_0_2_3.txt", header=TRUE)
raw7 <- read.csv("data/raw_data_turk7.csv")
lpsub <- subset(raw7, pic_name=="plot_large_p_small_n_30_100_0_2_3.png")
dt <- data.frame(xtabs(~response_no, data=lpsub))
names(dt) <- c(".sample", "picks")
picks <- data.frame(.sample=1:20)
picks <- merge(picks, dt, by=".sample", all.x=TRUE)
picks$picks[is.na(picks$picks)] <- 0

qplot(X1, X2, colour=factor(cl), facets=~.sample, data=lpt7, size=I(3)) + 
  theme_lineup +
  theme(legend.position="none") +
  geom_text(aes(label=picks), x=2.4, y=-1.8, data=picks, colour="grey70", 
            size=8, hjust=1)




# densities

ymax <- max(density(dists$dbn)$y)

ggplot() + 
  geom_density(data = dists, aes(x = dbn), 
               fill = "grey80", col = "grey80" ) +
      xlab("Binned (5,5) density") + ylab("")  +
  geom_segment(data = subset(dd3, len != 19), 
             aes(x = dbn.mean, xend = dbn.mean, 
                 y = rep(0.005*ymax,19), yend = rep(0.1*ymax,19)), 
             size=1, alpha = I(0.7)) + 
  geom_segment(data = subset(dd3, len == 19), 
             aes(x = dbn.mean, xend = dbn.mean, 
                 y = 0.005*ymax, yend = 0.2*ymax), 
             size=1, alpha = I(0.7), colour="darkorange") +
  geom_text(data = subset(dd3, len != 19), y = - 0.03*ymax, size = 2.5, 
            aes(x = dbn.mean, label = pos.1)) +
  geom_text(data = subset(dd3, len == 19), y = 0.25*ymax, size = 3, 
            colour="darkorange", aes(x = dbn.mean, label = pos.1)) 


ggsave("bin-dist-largep-5-5.pdf", height = 5, width = 5.5)        

### Sep dist

ymax <- max(density(dists$das)$y)

ggplot() + 
  geom_density(data = dists, aes(x = das), 
               fill = "grey80", col = "grey80" ) +
  xlab("Average separation density") + ylab("")  +
  geom_segment(data = subset(dd3, len != 19), 
               aes(x = das.mean, xend = das.mean, 
                   y = rep(0.005*ymax,19), yend = rep(0.1*ymax,19)), 
               size=1, alpha = I(0.7)) +
  geom_segment(data = subset(dd3, len == 19), 
               aes(x = das.mean, xend = das.mean, 
                   y = 0.005*ymax, yend = 0.2*ymax), 
               size=1, alpha = I(0.7), colour="darkorange") +
  geom_text(data = subset(dd3, len != 19), y = - 0.03*ymax, size = 2.5, 
            aes(x = das.mean, label = pos.1)) +
  geom_text(data = subset(dd3, len == 19), y = 0.25*ymax, size = 3, 
            colour="darkorange", aes(x = das.mean, label = pos.1)) 
  

# minimal separation
ymax <- max(density(dists$dms)$y)

ggplot() + 
  geom_density(data = dists, aes(x = dms), 
               fill = "grey80", col = "grey80" ) +
  xlab("Minimum Separation density") + ylab("")  +
  geom_segment(data = subset(dd3, len != 19), 
               aes(x = dms.mean, xend = dms.mean, 
                   y = rep(0.005*ymax,19), yend = rep(0.1*ymax,19)), 
               size=1, alpha = I(0.7)) +
  geom_segment(data = subset(dd3, len == 19), 
               aes(x = dms.mean, xend = dms.mean, 
                   y = 0.005*ymax, yend = 0.2*ymax), 
               size=1, alpha = I(0.7), colour="darkorange") +
  geom_text(data = subset(dd3, len != 19), y = - 0.03*ymax, size = 2.5, 
            aes(x = dms.mean, label = pos.1)) +
  geom_text(data = subset(dd3, len == 19), y = 0.25*ymax, size = 3, 
            colour="darkorange", aes(x = dms.mean, label = pos.1)) 

#####################
# minimal sep

ymax <- max(density(dists$dmin)$y)

ggplot() + 
  geom_density(data = dists, aes(x = dmin), 
               fill = "grey80", col = "grey80" ) +
  xlab("Minimum separation density") + ylab("")  +
  geom_segment(data = subset(dd3, len != 19), 
               aes(x = dmin.mean, xend = dmin.mean, 
                   y = rep(0.005*ymax,19), yend = rep(0.1*ymax,19)), 
               size=1, alpha = I(0.7)) + 
  geom_segment(data = subset(dd3, len == 19), 
               aes(x = dmin.mean, xend = dmin.mean, 
                   y = 0.005*ymax, yend = 0.2*ymax), 
               size=1, alpha = I(0.7), colour="darkorange") +
  geom_text(data = subset(dd3, len != 19), y = - 0.03*ymax, size = 2.5, 
            aes(x = dmin.mean, label = pos.1)) +
  geom_text(data = subset(dd3, len == 19), y = 0.25*ymax, size = 3, 
            colour="darkorange", aes(x = dmin.mean, label = pos.1)) 


##### 
# Dunn

ymax <- max(density(dists$ddunn)$y)

ggplot() + 
  geom_density(data = dists, aes(x = ddunn), 
               fill = "grey80", col = "grey80" ) +
  xlab("Dunn density") + ylab("")  +
  geom_segment(data = subset(dd3, len != 19), 
               aes(x = ddunn.mean, xend = ddunn.mean, 
                   y = rep(0.005*ymax,19), yend = rep(0.1*ymax,19)), 
               size=1, alpha = I(0.7)) + 
  geom_segment(data = subset(dd3, len == 19), 
               aes(x = ddunn.mean, xend = ddunn.mean, 
                   y = 0.005*ymax, yend = 0.2*ymax), 
               size=1, alpha = I(0.7), colour="darkorange") +
  geom_text(data = subset(dd3, len != 19), y = - 0.03*ymax, size = 2.5, 
            aes(x = ddunn.mean, label = pos.1)) +
  geom_text(data = subset(dd3, len == 19), y = 0.25*ymax, size = 3, 
            colour="darkorange", aes(x = ddunn.mean, label = pos.1)) 
