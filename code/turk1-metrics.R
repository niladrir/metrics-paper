
# re-calculating turk1-metrics
library(nullabor)

files <- dir("data/turk1")
files <- files[grep(".txt", files)]
results <- NULL

for(k in files){	
  ### Reading the lineup data
  dat <- read.table(sprintf("data/turk1/%s",k), header = T)
  
  ### Melting the data
  dat.m <- melt(dat, id = c("age", "grp", "weight"))
  
  
  ### Changing the categorical variable to a numerical variable
  dat.m$x <- as.numeric(dat.m$grp) 
  
  ### Breaking the variable name to get the position and type of plot (null or obs)
  dat.m$plot <- substring(dat.m$variable, 1, 3)
  dat.m$position <- substring(dat.m$variable, 4, 5)
  
  ### Finding the observed data
  obs <- dat.m[dat.m$plot == "obs", c("x", "value") ]
  
  ### Storing the pic name
  file.name <- k
  split.file.name <- unlist(strsplit(file.name,"\\."))
  split.2 <- unlist(strsplit(split.file.name[1],"\\_"))
  pic_name <- paste("plot_", split.2[2],"_", split.2[3], "_", split.2[4], "_", split.2[5], "_", split.2[6], ".png", sep = "" )
  
  ### Calculating the distance metrics on the lineups. Only binned, weighted bin, canberra distance 
  ### and hausdorff are calculated.
  
  dat.m <- dat.m[, c("x", "value", "position")]
  names(dat.m) <- c("group", "val", ".sample")
  lineup.dat <- dat.m
  lineup.dat$group <- as.factor(lineup.dat$group)

  metrics <- expand.grid(pos.1 = 1:20, pos.2 = 1:20)
  metrics$box.dist <- NA
  metrics$bin.dist <- NA
  
  for (i in 1:nrow(metrics)) {
    if (metrics$pos.1[i] == metrics$pos.2[i]) {
      metrics$box.dist[i] <- 0
      metrics$bin.dist[i] <- 0
    } else {
    .P <- subset(lineup.dat, .sample == metrics$pos.1[i])[,1:2]
    .PX <- subset(lineup.dat, .sample == metrics$pos.2[i])[,1:2]
    metrics$box.dist[i] <- box_dist(X=.P, PX=.PX)
    metrics$bin.dist[i] <- bin_dist(X=.P, PX=.PX, lineup.dat=lineup.dat,  X.bin = 2, Y.bin = 8)
    }
  }  
  res <- data.frame(pic_name, metrics)
  results <- rbind(results, res)
}

write.table(results, "turk1-metrics.txt", row.names = F)
