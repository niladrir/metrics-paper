files.txt <- dir("data/turk2","*.txt")

metrics <- NULL
for(k in 1:length(files.txt)){
  dat <- read.table(sprintf("data/turk2/%s",files.txt[k]), header = T)
  dat.m <- melt(dat, id.var=1)
  dat.m$.sample <- substring(dat.m$variable, 2)
  lineup.dat <- data.frame(x = dat.m$X, z = dat.m$value, .sample = dat.m$.sample)
  
  metrics.dat <- expand.grid(pos.1 = 1:20, pos.2 = 1:20)
  metrics.dat$bin.dist <- NA
  metrics.dat$reg.bin <- NA
  metrics.dat$reg.no.int <- NA
  
  for (i in 1:nrow(metrics.dat)) {
      .X <- subset(lineup.dat, .sample == metrics.dat$pos.1[i])[,1:2]
      .PX <- subset(lineup.dat, .sample == metrics.dat$pos.2[i])[,1:2]

      metrics.dat$bin.dist[i] <- bin_dist(X=.X, PX=.PX, lineup.dat=lineup.dat,  X.bin = 2, Y.bin = 2)
      metrics.dat$reg.bin[i] <- reg_dist(X=.X, PX=.PX) 
      metrics.dat$reg.no.int[i] <- reg_dist(X=.X, PX=.PX, intercept=FALSE) 
  }  
  
  
#  metrics.dat <- ddply(dat.pos, .(pos.1, pos.2), summarize, reg.bin = reg_bin_indx(pos.1, pos.2), bin.dist = bdist_mod_indx(pos.1, pos.2, nbin.X = 2, nbin.Y = 2), reg.no.int = reg_no_int_indx(pos.1, pos.2))
  pngname <- gsub("txt", "png", gsub("dat","plot", files.txt[k]))
  metrics.dat <- data.frame(metrics.dat, pic_name = pngname)
  metrics <- rbind(metrics, metrics.dat)
}

write.table(metrics, "turk2-metrics.txt", row.names = F)
